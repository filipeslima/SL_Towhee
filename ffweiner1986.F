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
      subroutine twh_ffweiner1986(lfailure,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes,nimprotypes,nhbondtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the weiner1986   *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 07-06-2005 by M.G. Martin                   *
c     * last modified 07-21-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
     &        ,nimprotypes,nhbondtypes
c     --- local variables
      integer ii,type,isame
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision rstar,epsilon
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)
      integer nislocal 
      character*(FFNAMELEN) ninlocal 
      dimension ninlocal(MAXSAMEIMP,MAXIMPNAME)
      character*(FFNAMELEN) hbnlocal
      dimension hbnlocal(MAXHBONDNAME)
      double precision hbclocal
      dimension hbclocal(MINHBONDCOEFF:MAXHBONDCOEFF)
      double precision onepi

      ffname = 'Weiner1986'
      call twh_classical_potential(GLB_SET,'12-6 plus 12-10 H-bond')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

      type = 0
c     ******************************************************************
c     * Weiner et al 1986 force field                                  *
c     ******************************************************************
c     --- an all-atom addition to the following united-atom force field
c     S.J. Weiner; P.A. Kollman; D.A. Case; U.C. Singh; C. Ghio;
c     G. Alagona; S. Profeta, Jr.; P. Weiner;
c     "A New Force Field for Molecular Mechanical Simulation of 
c     Nucleic Acids and Proteins", J. Am. Chem. Soc. 106, 765-784 (1984)
c     --- new all-atom types, plus some modifications to old parameters,
c     --- described in
c     S.J. Weiner; P.A. Kollman; D.T. Nguyen; D.A. Case;
c     "An All Atom Force Field for Simulations of Proteins and Nucleic 
c     Acids", J. Comp. Chem. 7, 230-252 (1986)

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C*'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.92d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 2.0d0 * 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 2.00d0
      epsilon = 0.15d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 3.0d0 * 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CB'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CD'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CE'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CF'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CG'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.09d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CI'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CJ'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- this one is new compared with weiner1984
      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CK'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CN'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0 + 1.0079d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CP'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- this one is new compared to weiner1984
      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.120d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CQ'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- this one is new compared to weiner1984
      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.120d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CR'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.80d0
      epsilon = 0.06d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- this one is new compared to weiner1984
      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.120d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CV'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- this one is new compared to weiner1984
      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.120d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CW'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.00d0
      epsilon = 0.02d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.00d0
      epsilon = 0.02d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.00d0
      epsilon = 0.02d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- this atom type remains from weiner1984, but the values changed
      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.54d0
      epsilon = 0.010d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.00d0
      epsilon = 0.02d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.00d0
      epsilon = 0.02d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.20d0
      epsilon = 0.016d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 0.0d0
      elelocal =  'LP'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'LP'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.75d0
      epsilon = 0.16d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.75d0
      epsilon = 0.16d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N*'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.75d0
      epsilon = 0.16d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.08d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.75d0
      epsilon = 0.16d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.75d0
      epsilon = 0.16d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NB'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.75d0
      epsilon = 0.16d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- new atom type compared to weiner1984
      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.75d0
      epsilon = 0.16d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NP'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.85d0
      epsilon = 0.12d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.60d0
      epsilon = 0.20d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.60d0
      epsilon = 0.20d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.65d0
      epsilon = 0.15d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 1.65d0
      epsilon = 0.15d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 2.10d0
      epsilon = 0.20d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 30.974d0
      elelocal =  ' P'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'P'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 2.00d0
      epsilon = 0.20d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 32.066d0
      elelocal =  ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      rstar = 2.00d0
      epsilon = 0.20d0
      nbclocal(1) = rstar*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = epsilon*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 32.066d0
      elelocal =  ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'SH'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return
      
c     ******************************************************************
c     * harmonic bonds                                                 *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-C2
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'C2'
c     C-C3
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'C3'
c     C-CH
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'CH'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.419d0
      vclocal(1) = (447.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-CB
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'CB'
c     CB-CN
      isame = isame + 1
      vnlocal(isame,1) = 'CB'
      vnlocal(isame,2) = 'CN'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.400d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-CD
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'CD'
c     CA-CD
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'CD'
c     CB-CD
      isame = isame + 1
      vnlocal(isame,1) = 'CB'
      vnlocal(isame,2) = 'CD'
c     CD-CD
      isame = isame + 1
      vnlocal(isame,1) = 'CD'
      vnlocal(isame,2) = 'CD'
c     CD-CN
      isame = isame + 1
      vnlocal(isame,1) = 'CD'
      vnlocal(isame,2) = 'CN'
c     C-CA [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'CA'
c     CA-CA [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'CA'
c     CA-CN [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'CN'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.444d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-CJ
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'CJ'
c     C-CM
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'CM'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (490.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-N
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'N'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.383d0
      vclocal(1) = (424.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-N*
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'N*'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.388d0
      vclocal(1) = (418.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-NA
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'NA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.358d0
      vclocal(1) = (457.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-NC
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'NC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.229d0
      vclocal(1) = (570.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-O
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'O'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.25d0
      vclocal(1) = (656.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-O2
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'O2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.364d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-OH
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'OH'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.495d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C2-C*
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'C*'
c     C*-CT [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'C*'
      vnlocal(isame,2) = 'CT'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C2-C2
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'C2'
c     C2-C3
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'C3'
c     CH-C2
      isame = isame + 1
      vnlocal(isame,1) = 'CH'
      vnlocal(isame,2) = 'C2'
c     CH-C3
      isame = isame + 1
      vnlocal(isame,1) = 'CH'
      vnlocal(isame,2) = 'C3'
c     CH-CH
      isame = isame + 1
      vnlocal(isame,1) = 'CH'
      vnlocal(isame,2) = 'CH'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.449d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C2-N
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'N'
c     C3-N
      isame = isame + 1
      vnlocal(isame,1) = 'C3'
      vnlocal(isame,2) = 'N'
c     CH-N
      isame = isame + 1
      vnlocal(isame,1) = 'CH'
      vnlocal(isame,2) = 'N'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.463d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C2-N2
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'N2'
c     C3-N2
      isame = isame + 1
      vnlocal(isame,1) = 'C3'
      vnlocal(isame,2) = 'N2'
c     CT-N2 [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'N2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.471d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C2-N3
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'N3'
c     C2-NT
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'NT'
c     C3-N3
      isame = isame + 1
      vnlocal(isame,1) = 'C3'
      vnlocal(isame,2) = 'N3'
c     CH-NT
      isame = isame + 1
      vnlocal(isame,1) = 'CH'
      vnlocal(isame,2) = 'NT'
c     CT-N3 [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'N3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.425d0
      vclocal(1) = (386.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C2-OH
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'OH'
c     CH-OH
      isame = isame + 1
      vnlocal(isame,1) = 'CH'
      vnlocal(isame,2) = 'OH'
c     C3-OH [new compared to weiner1984
      isame = isame + 1
      vnlocal(isame,1) = 'C3'
      vnlocal(isame,2) = 'OH'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.425d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C2-OS
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'OS'
c     C3-OS
      isame = isame + 1
      vnlocal(isame,1) = 'C3'
      vnlocal(isame,2) = 'OS'
c     CH-OS
      isame = isame + 1
      vnlocal(isame,1) = 'CH'
      vnlocal(isame,2) = 'OS'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C2-S
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'S'
c     C2-SH
      isame = isame + 1
      vnlocal(isame,1) = 'C2'
      vnlocal(isame,2) = 'SH'
c     C3-S
      isame = isame + 1
      vnlocal(isame,1) = 'C3'
      vnlocal(isame,2) = 'S'
c     CT-S [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'S'
c     CT-SH [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'SH'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.51d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C3-CM
      isame = isame + 1
      vnlocal(isame,1) = 'C3'
      vnlocal(isame,2) = 'CM'
c     CA-C2
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'C2'
c     CA-CT [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'CT'
c     CM-CT [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CM'
      vnlocal(isame,2) = 'CT'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.475d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C3-N*
      isame = isame + 1
      vnlocal(isame,1) = 'C3'
      vnlocal(isame,2) = 'N*'
c     CH-N*
      isame = isame + 1
      vnlocal(isame,1) = 'CH'
      vnlocal(isame,2) = 'N*'
c     CT-N*
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'N*'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.404d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CA-CB
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'CB'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.433d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CA-CJ
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'CJ'
c     CA-CM
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'CM'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.340d0
      vclocal(1) = (481.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CA-N2
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'N2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.381d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CA-NA
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'NA'
c     CG-NA
      isame = isame + 1
      vnlocal(isame,1) = 'CG'
      vnlocal(isame,2) = 'NA'
c     CW-NA [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CW'
      vnlocal(isame,2) = 'NA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.339d0
      vclocal(1) = (483.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CA-NC
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'NC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.459d0
      vclocal(1) = (388.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CB-C*
      isame = isame + 1
      vnlocal(isame,1) = 'CB'
      vnlocal(isame,2) = 'C*'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.370d0
      vclocal(1) = (520.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CB-CB
      isame = isame + 1
      vnlocal(isame,1) = 'CB'
      vnlocal(isame,2) = 'CB'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.374d0
      vclocal(1) = (436.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CB-N*
      isame = isame + 1
      vnlocal(isame,1) = 'CB'
      vnlocal(isame,2) = 'N*'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.391d0
      vclocal(1) = (414.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CB-NC
      isame = isame + 1
      vnlocal(isame,1) = 'CB'
      vnlocal(isame,2) = 'NB'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.354d0
      vclocal(1) = (461.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CB-NC
      isame = isame + 1
      vnlocal(isame,1) = 'CB'
      vnlocal(isame,2) = 'NC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.504d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CC-C2
      isame = isame + 1
      vnlocal(isame,1) = 'CC'
      vnlocal(isame,2) = 'C2'
c     CC-CT [new compared to weiner984]
      isame = isame + 1
      vnlocal(isame,1) = 'CC'
      vnlocal(isame,2) = 'CT'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.375d0
      vclocal(1) = (512.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CC-CF
      isame = isame + 1
      vnlocal(isame,1) = 'CC'
      vnlocal(isame,2) = 'CF'
c     CC-CV [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CC'
      vnlocal(isame,2) = 'CV'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (518.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CC-CG
      isame = isame + 1
      vnlocal(isame,1) = 'CC'
      vnlocal(isame,2) = 'CG'
c     CC-CW [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CC'
      vnlocal(isame,2) = 'CW'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.385d0
      vclocal(1) = (422.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CC-NA
      isame = isame + 1
      vnlocal(isame,1) = 'CC'
      vnlocal(isame,2) = 'NA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.394d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CC-NB
      isame = isame + 1
      vnlocal(isame,1) = 'CC'
      vnlocal(isame,2) = 'NB'
c     CF-NB
      isame = isame + 1
      vnlocal(isame,1) = 'CF'
      vnlocal(isame,2) = 'NB'
c     CV-NB
      isame = isame + 1
      vnlocal(isame,1) = 'CV'
      vnlocal(isame,2) = 'NB'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (440.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CE-N*
      isame = isame + 1
      vnlocal(isame,1) = 'CE'
      vnlocal(isame,2) = 'N*'
c     CK-N* [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CK'
      vnlocal(isame,2) = 'N*'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.304d0
      vclocal(1) = (529.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CE-NB
      isame = isame + 1
      vnlocal(isame,1) = 'CE'
      vnlocal(isame,2) = 'NB'
c     CK-NB [weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CK'
      vnlocal(isame,2) = 'NB'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.352d0
      vclocal(1) = (546.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CG-C*
      isame = isame + 1
      vnlocal(isame,1) = 'CG'
      vnlocal(isame,2) = 'C*'
c     C*-CW [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'C*'
      vnlocal(isame,2) = 'CW'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.324d0
      vclocal(1) = (502.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CI-NC
      isame = isame + 1
      vnlocal(isame,1) = 'CI'
      vnlocal(isame,2) = 'NC'
c     CQ-NC [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CQ'
      vnlocal(isame,2) = 'NC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.350d0
      vclocal(1) = (549.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CJ-CJ
      isame = isame + 1
      vnlocal(isame,1) = 'CJ'
      vnlocal(isame,2) = 'CJ'
c     CJ-CM [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CJ'
      vnlocal(isame,2) = 'CM'
c     CM-CM [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CM'
      vnlocal(isame,2) = 'CM'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.343d0
      vclocal(1) = (560.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CJ-CM
      isame = isame + 1
      vnlocal(isame,1) = 'CJ'
      vnlocal(isame,2) = 'CM'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (448.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CJ-N*
      isame = isame + 1
      vnlocal(isame,1) = 'CJ'
      vnlocal(isame,2) = 'N*'
c     CM-N* [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CM'
      vnlocal(isame,2) = 'N*'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.38d0
      vclocal(1) = (428.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CN-NA
      isame = isame + 1
      vnlocal(isame,1) = 'CN'
      vnlocal(isame,2) = 'NA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.343d0
      vclocal(1) = (477.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CP-NA
      isame = isame + 1
      vnlocal(isame,1) = 'CP'
      vnlocal(isame,2) = 'NA'
c     CR-NA [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CR'
      vnlocal(isame,2) = 'NA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (488.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CP-NB
      isame = isame + 1
      vnlocal(isame,1) = 'CP'
      vnlocal(isame,2) = 'NB'
c     CR-NB
      isame = isame + 1
      vnlocal(isame,1) = 'CR'
      vnlocal(isame,2) = 'NB'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (310.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CT-CT
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'CT'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.09d0
      vclocal(1) = (331.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CT-HC
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'HC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.01d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     H-N
      isame = isame + 1
      vnlocal(isame,1) = 'H'
      vnlocal(isame,2) = 'N'
c     H-N2
      isame = isame + 1
      vnlocal(isame,1) = 'H'
      vnlocal(isame,2) = 'N2'
c     H-NA
      isame = isame + 1
      vnlocal(isame,1) = 'H'
      vnlocal(isame,2) = 'NA'
c     H2-N
      isame = isame + 1
      vnlocal(isame,1) = 'H2'
      vnlocal(isame,2) = 'N'
c     H2-N2
      isame = isame + 1
      vnlocal(isame,1) = 'H2'
      vnlocal(isame,2) = 'N2'
c     H2-NT
      isame = isame + 1
      vnlocal(isame,1) = 'H2'
      vnlocal(isame,2) = 'NT'
c     H3-N2
      isame = isame + 1
      vnlocal(isame,1) = 'H3'
      vnlocal(isame,2) = 'N2'
c     H3-N3
      isame = isame + 1
      vnlocal(isame,1) = 'H3'
      vnlocal(isame,2) = 'N3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.96d0
      vclocal(1) = (553.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     HO-OH
      isame = isame + 1
      vnlocal(isame,1) = 'HO'
      vnlocal(isame,2) = 'OH'
c     HO-OS
      isame = isame + 1
      vnlocal(isame,1) = 'HO'
      vnlocal(isame,2) = 'OS'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.336d0
      vclocal(1) = (274.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     HS-SH
      isame = isame + 1
      vnlocal(isame,1) = 'HS'
      vnlocal(isame,2) = 'SH'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.679d0
      vclocal(1) = (600.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     LP-S
      isame = isame + 1
      vnlocal(isame,1) = 'LP'
      vnlocal(isame,2) = 'S'
c     LP-SH
      isame = isame + 1
      vnlocal(isame,1) = 'LP'
      vnlocal(isame,2) = 'SH'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.48d0
      vclocal(1) = (525.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     O2-P
      isame = isame + 1
      vnlocal(isame,1) = 'O2'
      vnlocal(isame,2) = 'P'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.61d0
      vclocal(1) = (230.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     OH-P
      isame = isame + 1
      vnlocal(isame,1) = 'OH'
      vnlocal(isame,2) = 'P'
c     OS-P
      isame = isame + 1
      vnlocal(isame,1) = 'OS'
      vnlocal(isame,2) = 'P'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.038d0
      vclocal(1) = (166.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     S-S
      isame = isame + 1
      vnlocal(isame,1) = 'S'
      vnlocal(isame,2) = 'S'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (335.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C-CT [bond remains from weiner1984, but values changed]
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'CT'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.080d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C*-HC [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'C*'
      vnlocal(isame,2) = 'HC'
c     CA-HC [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'HC'
c     CK-HC [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CK'
      vnlocal(isame,2) = 'HC'
c     CM-HC [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CM'
      vnlocal(isame,2) = 'HC'
c     CQ-HC [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CQ'
      vnlocal(isame,2) = 'HC'
c     CR-HC [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CR'
      vnlocal(isame,2) = 'HC'
c     CV-HC [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CV'
      vnlocal(isame,2) = 'HC'
c     CW-HC [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CW'
      vnlocal(isame,2) = 'HC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.449d0
      vclocal(1) = (355.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CT-N [bond remains from weiner1984, but values changed]
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'N'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.410d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CT-OH [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'OH'
c     CT-OS [new compared to weiner1984]
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'OS'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return
       
c     ******************************************************************
c     * harmonic angles                                                *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.7d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C*-CG-NA
      isame = isame + 1
      bnlocal(isame,1) = 'C*'
      bnlocal(isame,2) = 'CG'
      bnlocal(isame,3) = 'NA'
c     C*-CW-NA
      isame = isame + 1
      bnlocal(isame,1) = 'C*'
      bnlocal(isame,2) = 'CW'
      bnlocal(isame,3) = 'NA'
c     CG-CC-NA [very slightly changed from weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CG'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NA'
c     CW-CC-NA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CW'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.4d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-C2-C2
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'C2'
c     C-C2-CH
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'CH'
c     C2-C2-C2
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'C2'
c     CH-C2-C2
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'C2'
c     CH-C2-C3
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'C3'
c     CH-C2-CH
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'CH'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.3d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-C2-N
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'N'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.2d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-CB-CB
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'CB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-CB-NB
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'NB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.1d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-CH-C2
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'C2'
c     C-CH-C3
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'C3'
c     C-CH-CH
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'CH'
c     C-CT-CT
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.1d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-CH-N
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'N'
c     C-CT-N [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'N'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.7d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-CJ-CJ
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CJ'
      bnlocal(isame,3) = 'CJ'
c     C-CM-CJ
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'CJ'
c     C-CM-CM [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'CJ'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.7d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-CM-C3
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'C3'
c     CJ-CM-C3
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'C3'
c     CJ-CM-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.6d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-N*-CH
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CH'
c     C-N*-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.6d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-N*-CJ
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CJ'
c     C-N*-CM [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CM'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.9d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-N-C2
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'C2'
c     C-N-C3
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'C3'
c     C-N-CH
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'CH'
c     C-N-CT
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.8d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-N-H
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-N-H2
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'H2'
c     CA-N2-H
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'H'
c     CA-N2-H2
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'H2'
c     CA-N2-H3
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'H3'
c     H-N-H
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'H'
c     H2-N2-H2
      isame = isame + 1
      bnlocal(isame,1) = 'H2'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'H2'
c     H3-N2-H3
      isame = isame + 1
      bnlocal(isame,1) = 'H3'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'H3'
c     H3-N-H3
      isame = isame + 1
      bnlocal(isame,1) = 'H3'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'H3'
c     C-CA-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'HC'
c     C*-CW-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C*'
      bnlocal(isame,2) = 'CW'
      bnlocal(isame,3) = 'HC'
c     CA-CA-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'HC'
c     CB-CA-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'HC'
c     CC-CV-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'CV'
      bnlocal(isame,3) = 'HC'
c     CC-CW-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'CW'
      bnlocal(isame,3) = 'HC'
c     CN-CA-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CN'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'HC'
c     HC-CR-NA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CR'
      bnlocal(isame,3) = 'NA'
c     HC-CR-NB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CR'
      bnlocal(isame,3) = 'NB'
c     HC-CV-NB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CV'
      bnlocal(isame,3) = 'NB'
c     HC-CW-NA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CW'
      bnlocal(isame,3) = 'NA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.4d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-NA-C
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.2d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-NA-CA
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.8d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-NA-H
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.5d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-NC-CA
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'NC'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.0d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-OH-HO
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'OH'
      bnlocal(isame,3) = 'HO'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.6d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-C-N
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'N'
c     C3-C-N
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'N'
c     CH-C-N
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'N'
c     CT-C-N
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'N'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.4d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     C3-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     CH-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     CT-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-C-O2
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O2'
c     C3-C-O2
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O2'
c     CT-C-O2 [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.2d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-C2-N
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'N'
c     C2-C2-N2
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'N2'
c     C2-C2-N3
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'N3'
c     NT-C2-C
      isame = isame + 1
      bnlocal(isame,1) = 'NT'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'C'
c     NT-C2-C2
      isame = isame + 1
      bnlocal(isame,1) = 'NT'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'C2'
c     CT-CT-N2 [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'N2'
c     CT-CT-N3 [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'N3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-C2-OS
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'OS'
c     C2-CH-N*
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'N*'
c     C2-CH-OH
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'OH'
c     C2-CH-OS
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'OS'
c     C3-C2-OS
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'OS'
c     C3-CH-N
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'N'
c     C3-CH-OH
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'OH'
c     CH-C2-OH
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'OH'
c     CH-C2-OS
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'OS'
c     CH-CH-N*
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'N*'
c     CH-CH-OH
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'OH'
c     CH-CH-OS
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'OS'
c     N*-CH-OS
      isame = isame + 1
      bnlocal(isame,1) = 'N*'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'OS'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.7d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-C2-S
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'S'
c     CH-C2-S
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'S'
c     CT-CT-S
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'S'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 131.9d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-CC-CF
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CF'
c     C2-CC-CV [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CV'
c     CT-CC-CV [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CV'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-CC-CG [modified very slightly from weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CG'
c     C2-CC-CW [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CW'
c     CT-CC-CW [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CW'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.2d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-CC-NA
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NA'
c     CT-CC-NA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-CC-NB [very slightly changed from weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NB'
c     CT-CC-NB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.5d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-CH-C3
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'C3'
c     C3-CH-C3
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'C3'
c     CH-CH-C2
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'C2'
c     CH-CH-C3
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'C3'
c     CH-CH-CH
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'CH'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.4d0
      bclocal(1) = (38.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-N-H
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'H'
c     C3-N-H
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'H'
c     CH-N-H
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'H'
c     CT-N-H
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.4d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-N2-H2
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'H2'
c     C2-N2-H3
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'H3'
c     C3-N2-H2
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'H2'
c     CT-N2-H3 [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'H3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-N3-H3
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'N3'
      bnlocal(isame,3) = 'H3'
c     C2-NT-H2
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'NT'
      bnlocal(isame,3) = 'H2'
c     C3-N3-H3
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'N3'
      bnlocal(isame,3) = 'H3'
c     CH-NT-H2
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'NT'
      bnlocal(isame,3) = 'H2'
c     H2-NT-H2
      isame = isame + 1
      bnlocal(isame,1) = 'H2'
      bnlocal(isame,2) = 'NT'
      bnlocal(isame,3) = 'H2'
c     H3-N3-H3
      isame = isame + 1
      bnlocal(isame,1) = 'H3'
      bnlocal(isame,2) = 'N3'
      bnlocal(isame,3) = 'H3'
c     C-CT-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'HC'
c     C*-CT-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C*'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'HC'
c     CA-CT-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'HC'
c     CC-CT-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'HC'
c     CM-CT-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'HC'
c     CT-CT-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'HC'
c     CT-N3-H3 [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'N3'
      bnlocal(isame,3) = 'H3'
c     HC-CT-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'HC'
c     HC-CT-N* [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'N*'
c     HC-CT-N2 [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'N2'
c     HC-CT-N3 [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'NE'
c     HC-CT-OH [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'OH'
c     HC-CT-OS [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'OS'
c     HC-CT-S [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'S'
c     HC-CT-SH [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'SH'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.5d0
      bclocal(1) = (55.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-OH-HO
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'OH'
      bnlocal(isame,3) = 'HO'
c     C3-OH-HO
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'OH'
      bnlocal(isame,3) = 'HO'
c     C2-OS-HO
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'HO'
c     CH-OH-HO
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'OH'
      bnlocal(isame,3) = 'HO'
c     CH-OS-HO
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'HO'
c     CT-OH-HO [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'OH'
      bnlocal(isame,3) = 'HO'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.8d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-OS-C2
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'C2'
c     C2-OS-C3
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'C3'
c     CH-OS-CH
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'CH'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.5d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-OS-P
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'P'
c     C3-OS-P
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'P'
c     CH-OS-P
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'P'
c     CT-OS-P [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'P'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 98.9d0
      bclocal(1) = (62.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-S-C3
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'C3'
c     CT-S-CT [new compared to weiner1986]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.7d0
      bclocal(1) = (68.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-S-S
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'S'
c     C3-S-S
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'S'
c     CT-S-S [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'S'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 96.0d0
      bclocal(1) = (44.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C2-SH-HS
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'SH'
      bnlocal(isame,3) = 'HS'
c     C3-SH-HS
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'SH'
      bnlocal(isame,3) = 'HS'
c     CT-SH-HS [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'SH'
      bnlocal(isame,3) = 'HS'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.3d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CA-CB-CB
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'CB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 132.4d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CA-CB-NB
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'NB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CA-CD-CD
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CD'
c     CB-CD-CD
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CD'
c     CD-C-CD
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CD'
c     CD-CA-CD
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CD'
c     CD-CD-C
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'C'
c     CD-CD-CD
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CD'
c     CD-CD-CN
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CN'
c     C-CA-CA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     CA-C-CA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CA'
c     CA-CA-CA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     CA-CA-CB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CB'
c     CA-CA-CN [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CN'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.0d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CA-CJ-CJ
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CJ'
      bnlocal(isame,3) = 'CJ'
c     CA-CM-CM [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'CM'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.2d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CA-N2-C2
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'C2'
c     CA-N2-C3
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'C3'
c     CA-N2-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'N2'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.0d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CA-NA-H
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.2d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CA-NC-CB
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'NC'
      bnlocal(isame,3) = 'CB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.6d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CA-NC-CI
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'NC'
      bnlocal(isame,3) = 'CI'
c     NC-C-N*
      isame = isame + 1
      bnlocal(isame,1) = 'NC'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'N*'
c     CA-NC-CQ [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'NC'
      bnlocal(isame,3) = 'CQ'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.6d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-C*-C2
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'C2'
c     CB-C*-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.4d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-C*-CG
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'CG'
c     CB-C*-CW [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'CW'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.3d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-C-NA
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'NA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.8d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.5d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-CA-N2
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'N2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.3d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-CA-NC
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'NC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.2d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-CB-N*
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'N*'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.4d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-CB-NB
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'NB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.7d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-CB-NC
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'NC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.7d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-CN-CD
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'CN'
      bnlocal(isame,3) = 'CD'
c     CA-CN-CB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CN'
      bnlocal(isame,3) = 'CB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.4d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-CN-NA
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'CN'
      bnlocal(isame,3) = 'NA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.8d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-N*-C3
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'C3'
c     CB-N*-CH
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CH'
c     CB-N*-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.4d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-N*-CE
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CE'
c     CB-N*-CK [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CK'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.8d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-NB-CE
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'NB'
      bnlocal(isame,3) = 'CE'
c     CB-NB-CK [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'NB'
      bnlocal(isame,3) = 'CK'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-NC-CI
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'NC'
      bnlocal(isame,3) = 'CI'
c     CB-NC-CQ [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'NC'
      bnlocal(isame,3) = 'CQ'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.1d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CC-C2-CH
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'CH'
c     CC-CT-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.9d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CC-CF-NB
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'CF'
      bnlocal(isame,3) = 'NB'
c     CG-CC-NB
      isame = isame + 1
      bnlocal(isame,1) = 'CG'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NB'
c     CC-CV-NB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'CV'
      bnlocal(isame,3) = 'NB'
c     CW-CC-NB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CW'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.9d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CC-CG-NA
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'CG'
      bnlocal(isame,3) = 'NA'
c     CF-CC-NA
      isame = isame + 1
      bnlocal(isame,1) = 'CF'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NA'
c     CC-CW-NA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'CW'
      bnlocal(isame,3) = 'NA'
c     CV-CC-NA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CV'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.3d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CC-NA-CP
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'CP'
c     CP-NA-CG
      isame = isame + 1
      bnlocal(isame,1) = 'CP'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'CG'
c     CC-NA-CR [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'CR'
c     CG-NA-CR [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CG'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'CR'
c     CR-NA-CW [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CR'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'CW'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.3d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CC-NA-H [changed very slightly from weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'H'
c     CG-NA-H
      isame = isame + 1
      bnlocal(isame,1) = 'CG'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'H'
c     CP-NA-H [changed very slightly from weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CP'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'H'
c     CR-NA-H [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CR'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.3d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CC-NB-CP
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'NB'
      bnlocal(isame,3) = 'CP'
c     CF-NB-CP
      isame = isame + 1
      bnlocal(isame,1) = 'CF'
      bnlocal(isame,2) = 'NB'
      bnlocal(isame,3) = 'CP'
c     CC-NB-CR [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CC'
      bnlocal(isame,2) = 'NB'
      bnlocal(isame,3) = 'CR'
c     CF-NB-CR [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CF'
      bnlocal(isame,2) = 'NB'
      bnlocal(isame,3) = 'CR'
c     CR-NB-CV [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CR'
      bnlocal(isame,2) = 'NB'
      bnlocal(isame,3) = 'CV'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CD-C-OH
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'OH'
c     CD-CA-C2
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'C2'
c     N2-CA-N2
      isame = isame + 1
      bnlocal(isame,1) = 'N2'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'N2'
c     C2-CA-CA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     CA-C-OH [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'OH'
c     CA-CA-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 134.9d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CD-CB-C*
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'C*'
c     C*-CB-CA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C*'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.2d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CD-CB-CN
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'CN'
c     CA-CB-CN [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'CN'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 132.8d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CD-CN-NA
      isame = isame + 1
      bnlocal(isame,1) = 'CD'
      bnlocal(isame,2) = 'CN'
      bnlocal(isame,3) = 'NA'
c     CA-CN-NA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CN'
      bnlocal(isame,3) = 'NA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.8d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CE-N*-C3
      isame = isame + 1
      bnlocal(isame,1) = 'CE'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'C3'
c     CE-N*-CH
      isame = isame + 1
      bnlocal(isame,1) = 'CE'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CH'
c     C3-N*-CK [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CK'
c     CE-N*-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CE'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CT'
c     CH-N*-CK [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CK'
c     CK-N*-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CK'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CG-C*-C2
      isame = isame + 1
      bnlocal(isame,1) = 'CG'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'C2'
c     C2-C*-CW [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C2'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'CW'
c     CT-C*-CW [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'CW'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.0d0
      bclocal(1) = (65.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CH-C-O2
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CH-C-OH
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'OH'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.6d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CH-C2-C*
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'C*'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.0d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CH-C2-CA
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'CA'
c     CA-CT-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.6d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CH-C2-SH
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'C2'
      bnlocal(isame,3) = 'SH'
c     CT-CT-SH [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'SH'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.0d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CH-N-C2
      isame = isame + 1
      bnlocal(isame,1) = 'CH'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'C2'
c     CT-N-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.1d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CJ-C-NA
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'NA'
c     CM-C-NA
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'NA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.3d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CJ-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     CM-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.1d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CJ-CA-N2
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'N2'
c     CM-CA-N2 [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'N2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.5d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CJ-CA-NC
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'NC'
c     CM-CA-NC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'NC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.2d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CJ-CJ-N*
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'CJ'
      bnlocal(isame,3) = 'N*'
c     CJ-N*-CH
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CH'
c     CM-CJ-N*
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'CJ'
      bnlocal(isame,3) = 'N*'
c     CJ-N*-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CT'
c     CM-CM-N* [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'N*'
c     CM-N*-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.8d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CN-CB-C*
      isame = isame + 1
      bnlocal(isame,1) = 'CN'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'C*'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.6d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CN-NA-CG
      isame = isame + 1
      bnlocal(isame,1) = 'CN'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'CG'
c     NB-CP-NA
      isame = isame + 1
      bnlocal(isame,1) = 'NB'
      bnlocal(isame,2) = 'CP'
      bnlocal(isame,3) = 'NA'
c     CN-NA-CW [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CN'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'CW'
c     NA-CR-NB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'NA'
      bnlocal(isame,2) = 'CR'
      bnlocal(isame,3) = 'NB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.1d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CN-NA-H [slight change from weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CN'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.5d0
      bclocal(1) = (47.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     HO-OH-HO
      isame = isame + 1
      bnlocal(isame,1) = 'HO'
      bnlocal(isame,2) = 'OH'
      bnlocal(isame,3) = 'HO'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.5d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     HO-OH-P
      isame = isame + 1
      bnlocal(isame,1) = 'HO'
      bnlocal(isame,2) = 'OH'
      bnlocal(isame,3) = 'P'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 96.7d0
      bclocal(1) = (600.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     LP-S-C2
      isame = isame + 1
      bnlocal(isame,1) = 'LP'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'C2'
c     LP-S-C3
      isame = isame + 1
      bnlocal(isame,1) = 'LP'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'C3'
c     LP-S-S
      isame = isame + 1
      bnlocal(isame,1) = 'LP'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'S'
c     LP-SH-C2
      isame = isame + 1
      bnlocal(isame,1) = 'LP'
      bnlocal(isame,2) = 'SH'
      bnlocal(isame,3) = 'C2'
c     LP-SH-HS
      isame = isame + 1
      bnlocal(isame,1) = 'LP'
      bnlocal(isame,2) = 'SH'
      bnlocal(isame,3) = 'H2'
c     C3-SH-LP [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'SH'
      bnlocal(isame,3) = 'LP'
c     CT-S-LP [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'LP'
c     CT-SH-LP [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'SH'
      bnlocal(isame,3) = 'LP'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 160.0d0
      bclocal(1) = (600.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     LP-S-LP
      isame = isame + 1
      bnlocal(isame,1) = 'LP'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'LP'
c     LP-SH-LP
      isame = isame + 1
      bnlocal(isame,1) = 'LP'
      bnlocal(isame,2) = 'SH'
      bnlocal(isame,3) = 'LP'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.9d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     N*-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'N*'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.9d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     N*-CE-NB
      isame = isame + 1
      bnlocal(isame,1) = 'N*'
      bnlocal(isame,2) = 'CE'
      bnlocal(isame,3) = 'NB'
c     N*-CK-NB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'N*'
      bnlocal(isame,2) = 'CK'
      bnlocal(isame,3) = 'NB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.9d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     N-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.7d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     N-CH-C2
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'C2'
c     N-CH-CH
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'CH'
c     NT-CH-C
      isame = isame + 1
      bnlocal(isame,1) = 'NT'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'C'
c     NT-CH-C2
      isame = isame + 1
      bnlocal(isame,1) = 'NT'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'C2'
c     NT-CH-CH
      isame = isame + 1
      bnlocal(isame,1) = 'NT'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'CH'
c     C3-CH-NT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C3'
      bnlocal(isame,2) = 'CH'
      bnlocal(isame,3) = 'NT'
c     CT-CT-N [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'N'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.4d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     NA-C-N*
      isame = isame + 1
      bnlocal(isame,1) = 'NA'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'N*'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.6d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     NA-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'NA'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     NA-CA-N2
      isame = isame + 1
      bnlocal(isame,1) = 'NA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'N2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.3d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     NA-CA-NC
      isame = isame + 1
      bnlocal(isame,1) = 'NA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'NC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.7d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     NA-CP-NA [changed slightly from weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'NA'
      bnlocal(isame,2) = 'CP'
      bnlocal(isame,3) = 'NA'
c     NA-CR-NA [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'NA'
      bnlocal(isame,2) = 'CR'
      bnlocal(isame,3) = 'NA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.5d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     NC-C-O
      isame = isame + 1
      bnlocal(isame,1) = 'NC'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.3d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     NC-CA-N2 [changed slightly from weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'NC'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'N2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.0d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     NC-CB-N*
      isame = isame + 1
      bnlocal(isame,1) = 'NC'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'N*'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.1d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     NC-CI-NC
      isame = isame + 1
      bnlocal(isame,1) = 'NC'
      bnlocal(isame,2) = 'CI'
      bnlocal(isame,3) = 'NC'
c     NC-CQ-NC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'NC'
      bnlocal(isame,2) = 'CQ'
      bnlocal(isame,3) = 'NC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.0d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     O2-C-O2
      isame = isame + 1
      bnlocal(isame,1) = 'O2'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'O2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.9d0
      bclocal(1) = (140.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     O2-P-O2
      isame = isame + 1
      bnlocal(isame,1) = 'O2'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'O2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.2d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     OH-P-O2
      isame = isame + 1
      bnlocal(isame,1) = 'OH'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'O2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.2d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     OS-P-O2
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'O2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 102.6d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     OS-P-OH
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'OH'
c     OS-P-OS
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'OS'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.7d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-CM-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'CT'
c     CM-CM-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.7d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-CM-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'HC'
c     CM-CM-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'HC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.2d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-N*-H [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'H'
c     CJ-N*-H [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CJ'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'H'
c     CM-N*-H [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CM'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.3d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CA-CM-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'HC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.6d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-C*-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.8d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-C*-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'HC'
c     CW-C*-HC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CW'
      bnlocal(isame,2) = 'C*'
      bnlocal(isame,3) = 'HC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.3d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CB-N*-H [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CB'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'H'
c     CE-N*-H [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CE'
      bnlocal(isame,2) = 'N*'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (40.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CT-CT-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CT-CT-N* [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'N*'
c     CT-CT-OH [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'OH'
c     CT-CT-OS [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'OS'
c     N*-CT-OS [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'N*'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'OS'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (60.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CT-OS-CT [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CT'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.3d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CW-NA-H [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'CW'
      bnlocal(isame,2) = 'NA'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.0d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     HC-CK-N* [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CK'
      bnlocal(isame,3) = 'N*'
c     HC-CK-NB [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CK'
      bnlocal(isame,3) = 'NB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.10d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     HC-CM-N* [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'N*'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.4d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     HC-CQ-NC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CQ'
      bnlocal(isame,3) = 'NC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (38.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     HC-CT-N [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'N'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 92.10d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     HS-SH-HS [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'HS'
      bnlocal(isame,2) = 'SH'
      bnlocal(isame,3) = 'HS'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.2d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     N*-CB-NC [new compared to weiner1984]
      isame = isame + 1
      bnlocal(isame,1) = 'N*'
      bnlocal(isame,2) = 'CB'
      bnlocal(isame,3) = 'NC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Amber/Charmm/Gromos style torsions                             *
c     * 1-4 scale factor of 1/2 is mentioned in the Cornell Amber paper*
c     * and is used here                                               *
c     * The listed values are divided by the total number of torsions  *
c     * this is done in vtorsion
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-C2-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-CA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.4d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-CB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-CD-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CD'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-CH-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.1d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-CJ-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CJ'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.1d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-CM-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-CT-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (10.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-N-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.8d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-N*-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N*'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.4d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-NA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (8.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-NC-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NC'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.8d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C-OH-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'OH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C*-C2-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C*'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.8d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C*-CB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C*'
      tnlocal(isame,3) = 'CB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (23.6d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C*-CG-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C*'
      tnlocal(isame,3) = 'CG'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C*-CT-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C*'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (23.6d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C*-CW-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C*'
      tnlocal(isame,3) = 'CW'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-C2-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-CA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-CC-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-CH-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-N-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-N2-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'N2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.4d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-N3-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'N3'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-NT-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'NT'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.5d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-OH-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'OH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.45d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-OS-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-S-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.75d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-SH-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'SH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-CA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (10.2d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-CB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-CD-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CD'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.7d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-CJ-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CJ'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.7d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-CM-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (10.6d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-CN-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CN'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-CT-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.8d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-N2-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'N2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-NA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'NA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (9.6d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CA-NC-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'NC'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (16.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CB-CB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CB'
      tnlocal(isame,3) = 'CB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (20.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CB-CN-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CB'
      tnlocal(isame,3) = 'CN'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.6d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CB-N*-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CB'
      tnlocal(isame,3) = 'N*'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.1d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CB-NB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CB'
      tnlocal(isame,3) = 'NB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (14.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CC-CF-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CF'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (15.9d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CC-CG-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CG'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CC-CT-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (14.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CC-CV-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CV'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (15.9d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CC-CW-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CW'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.6d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CC-NA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'NA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.8d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CC-NB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'NB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CD-CD-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CD'
      tnlocal(isame,3) = 'CD'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CD-CN-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CD'
      tnlocal(isame,3) = 'CN'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.7d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CE-N*-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CE'
      tnlocal(isame,3) = 'N*'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (20.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CE-NB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CE'
      tnlocal(isame,3) = 'NB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.8d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CF-NB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CF'
      tnlocal(isame,3) = 'NB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CG-NA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CG'
      tnlocal(isame,3) = 'NA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CH-CH-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CH-N-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CH-N*-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'N*'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CH-NT-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'NT'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.5d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CH-OH-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'OH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.45d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CH-OS-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (13.5d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CI-NC-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CI'
      tnlocal(isame,3) = 'NC'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (24.4d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CJ-CJ-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CJ'
      tnlocal(isame,3) = 'CJ'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (24.4d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CJ-CM-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CJ'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (7.4d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CJ-N*-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CJ'
      tnlocal(isame,3) = 'N*'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.7d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CK-N*-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CK'
      tnlocal(isame,3) = 'N*'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (20.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CK-NB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CK'
      tnlocal(isame,3) = 'NB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (24.4d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CM-CM-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CM-CT-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (7.4d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CM-N*-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'N*'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (12.2d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CN-NA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN'
      tnlocal(isame,3) = 'NA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (9.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CP-NA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CP'
      tnlocal(isame,3) = 'NA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (10.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CP-NB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CP'
      tnlocal(isame,3) = 'NB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (13.5d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CQ-NC-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CQ'
      tnlocal(isame,3) = 'NC'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (9.3d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CR-NA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CR'
      tnlocal(isame,3) = 'NA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (10.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CR-NB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CR'
      tnlocal(isame,3) = 'NB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.3d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CT-CT-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CT-N-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CT-N*-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'N*'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CT-N2-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'N2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.4d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CT-N3-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'N3'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.5d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CT-OH-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'OH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.15d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CT-OS-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CT-S-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.75d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CT-SH-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'SH'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.8d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CV-NB-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CV'
      tnlocal(isame,3) = 'NB'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.0d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-CW-NA-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CW'
      tnlocal(isame,3) = 'NA'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.75d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-OH-P-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'OH'
      tnlocal(isame,3) = 'P'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.75d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X-OS-P-X
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'P'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.200d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     O-C-C2-N
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'N'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     O-C-CH-C2
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'C2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     O-C-CH-N
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'N'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     O-C-CH-CH
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'CH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (2.000d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-C2-C2-OH
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'OH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (2.000d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OH-C2-C2-OH
      isame = isame + 1
      tnlocal(isame,1) = 'OH'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'OH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (2.000d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-C2-C2-OS
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'OS'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (1.000d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-C2-CH-OS
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'OS'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (1.000d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-C2-CH-OH
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'OH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (1.000d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OH-C2-CH-OH
      isame = isame + 1
      tnlocal(isame,1) = 'OH'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'OH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C2-C2-S-LP
      isame = isame + 1
      tnlocal(isame,1) = 'C2'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'LP'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CH-C2-SH-LP
      isame = isame + 1
      tnlocal(isame,1) = 'CH'
      tnlocal(isame,2) = 'C2'
      tnlocal(isame,3) = 'SH'
      tnlocal(isame,4) = 'LP'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(3) = (1.000d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1) = (1.000d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-CH-C2-OH
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'OH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OH-CH-CH-OH
      isame = isame + 1
      tnlocal(isame,1) = 'OH'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'OH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-CH-CH-OH
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'OH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-CH-CH-OS
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'OS'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.710d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     HC-CM-CM-CT
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'CT'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.590d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C-CM-CM-HC
      isame = isame + 1
      tnlocal(isame,1) = 'C'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'HC'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.590d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     N*-CM-CM-CT
      isame = isame + 1
      tnlocal(isame,1) = 'N*'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'CT'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.590d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CA-CM-CM-HC
      isame = isame + 1
      tnlocal(isame,1) = 'CA'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'HC'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (9.510d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     N*-CM-CM-CA
      isame = isame + 1
      tnlocal(isame,1) = 'N*'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'CA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.710d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     HC-CM-CM-HC
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'HC'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (9.510d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     N*-CM-CM-C
      isame = isame + 1
      tnlocal(isame,1) = 'N*'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'C'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.590d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     N*-CM-CM-HC
      isame = isame + 1
      tnlocal(isame,1) = 'N*'
      tnlocal(isame,2) = 'CM'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'HC'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.067d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     N-CT-C-O
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'O'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.067d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     HC-CT-C-O
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'O'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.067d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CT-CT-C-O
      isame = isame + 1
      tnlocal(isame,1) = 'CT'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'O'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.200d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.383d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CT-OS-CT-CT
      isame = isame + 1
      tnlocal(isame,1) = 'CT'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'CT'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.144d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-CT-CT-OS
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'OS'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.144d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-CT-CT-OH
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'OH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.144d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OH-CT-CT-OH
      isame = isame + 1
      tnlocal(isame,1) = 'OH'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'CT'
      tnlocal(isame,4) = 'OH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.650d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 1.0d0
      tclocal(1) = (2.500d0)*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     H-N-C-O
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'N'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'O'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.725d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C2-OS-C2-C3
      isame = isame + 1
      tnlocal(isame,1) = 'C2'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'C3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (1.450d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C2-OS-C2-C2
      isame = isame + 1
      tnlocal(isame,1) = 'C2'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'C2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (1.450d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C3-OS-C2-C3
      isame = isame + 1
      tnlocal(isame,1) = 'C3'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'C2'
      tnlocal(isame,4) = 'C3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.725d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.100d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CH-OS-CH-C2
      isame = isame + 1
      tnlocal(isame,1) = 'CH'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'C2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.725d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CH-OS-CH-CH
      isame = isame + 1
      tnlocal(isame,1) = 'CH'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'CH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.725d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C2-OS-CH-C2
      isame = isame + 1
      tnlocal(isame,1) = 'C2'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'C2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.725d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.100d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C3-OS-CH-C3
      isame = isame + 1
      tnlocal(isame,1) = 'C3'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'C3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.725d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.000d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CH-OS-CH-N*
      isame = isame + 1
      tnlocal(isame,1) = 'CH'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'N*'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.725d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C2-OS-CH-C3
      isame = isame + 1
      tnlocal(isame,1) = 'C2'
      tnlocal(isame,2) = 'OS'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'C3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.750d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.250d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OH-P-OS-C3
      isame = isame + 1
      tnlocal(isame,1) = 'OH'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'C3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.250d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.750d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-P-OS-C2
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'C2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.750d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.250d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OH-P-OS-C2
      isame = isame + 1
      tnlocal(isame,1) = 'OH'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'C2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.250d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.750d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OH-P-OS-CT
      isame = isame + 1
      tnlocal(isame,1) = 'OH'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'CT'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.750d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.250d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-P-OS-CH
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'CH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.250d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.750d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OH-P-OS-CH
      isame = isame + 1
      tnlocal(isame,1) = 'OH'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'CH'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.250d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (0.750d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OH-P-OS-CT
      isame = isame + 1
      tnlocal(isame,1) = 'OH'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'CT'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.750d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(1+3) = (0.250d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     OS-P-OS-C3
      isame = isame + 1
      tnlocal(isame,1) = 'OS'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'C3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     LP-S-S-LP
      isame = isame + 1
      tnlocal(isame,1) = 'LP'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'LP'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     LP-S-S-C2
      isame = isame + 1
      tnlocal(isame,1) = 'LP'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'C2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.600d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (3.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C2-S-S-C2
      isame = isame + 1
      tnlocal(isame,1) = 'C2'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'C2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.600d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(1+3) = (3.500d0)*(CON_KCALTOK)
      tclocal(3+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2+3) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CT-S-S-CT
      isame = isame + 1
      tnlocal(isame,1) = 'CT'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'CT'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,19)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     LP-S-S-CT
      isame = isame + 1
      tnlocal(isame,1) = 'LP'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'CT'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * impcoeff(1)*[ 1 + Cos( impcoeff(2)* theta - impcoeff(3)mm) ]   *
c     ******************************************************************
c     --- They are listed as i-j-k-l where the torsion is the angle 
c     --- between i-k-l and j-k-l.  i,j, and l are bonded to k and not 
c     --- to each other
c     --- This torsion is only applied once per planar center even in 
c     --- the case of C-N-H2 where you could apply it twice.
c     --- when this is defined in towhee_input it is read in from 
c     --- atom k where ijimp2 = l: ijimp3 = i, ijimp4 = j.

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,4)
      iclocal(1) = (10.5d0)*(CON_KCALTOK)
      iclocal(3) = (180.0d0)*(onepi/180.0d0)
      iclocal(2) = 2.0d0
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     X-X-C-O
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'C'
      ninlocal(isame,4) = 'O'
c     X-N2-CA-N2
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'N2'
      ninlocal(isame,3) = 'CA'
      ninlocal(isame,4) = 'N2'
c     X-O2-C-O2
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'O2'
      ninlocal(isame,3) = 'C'
      ninlocal(isame,4) = 'O2'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,4)
      iclocal(1) = (1.0d0)*(CON_KCALTOK)
      iclocal(3) = (180.0d0)*(onepi/180.0d0)
      iclocal(2) = 2.0d0
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     X-X-N-H
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'N'
      ninlocal(isame,4) = 'H'
c     X-X-NA-H
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'NA'
      ninlocal(isame,4) = 'H'
c     X-CH-N-C2
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'CH'
      ninlocal(isame,3) = 'N'
      ninlocal(isame,4) = 'C2'
c     X-H2-N-H2
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'H2'
      ninlocal(isame,3) = 'N'
      ninlocal(isame,4) = 'H2'
c     X-CT-N-CT [new compared to weiner1984]
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'CT'
      ninlocal(isame,3) = 'N'
      ninlocal(isame,4) = 'CT'
c     X-X-N2-H3 [new compared to weiner1984]
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'N2'
      ninlocal(isame,4) = 'H3'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,4)
      iclocal(1) = (14.0d0)*(CON_KCALTOK)
      iclocal(3) = (180.0d0)*(onepi/180.0d0)
      iclocal(2) = 3.0d0
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     X-C2-CH-X
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'C2'
      ninlocal(isame,3) = 'CH'
      ninlocal(isame,4) = 'wild'
c     X-CH-CH-X
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'CH'
      ninlocal(isame,3) = 'CH'
      ninlocal(isame,4) = 'wild'
c     X-CH-N-C
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'CH'
      ninlocal(isame,3) = 'N'
      ninlocal(isame,4) = 'C'
c     C3-CH-NT-C [new from weiner1984]
      isame = isame + 1
      ninlocal(isame,1) = 'C3'
      ninlocal(isame,2) = 'CH'
      ninlocal(isame,3) = 'NT'
      ninlocal(isame,4) = 'C'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,4)
      iclocal(1) = (7.0d0)*(CON_KCALTOK)
      iclocal(3) = (180.0d0)*(onepi/180.0d0)
      iclocal(2) = 3.0d0
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     C2-CH-C-N3
      isame = isame + 1
      ninlocal(isame,1) = 'C2'
      ninlocal(isame,2) = 'CH'
      ninlocal(isame,3) = 'C'
      ninlocal(isame,4) = 'N3'
c     C3-CH-CA-C3
      isame = isame + 1
      ninlocal(isame,1) = 'C2'
      ninlocal(isame,2) = 'CH'
      ninlocal(isame,3) = 'CA'
      ninlocal(isame,4) = 'C3'
c     CH-CH-C-N3
      isame = isame + 1
      ninlocal(isame,1) = 'CH'
      ninlocal(isame,2) = 'CH'
      ninlocal(isame,3) = 'C'
      ninlocal(isame,4) = 'N3'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,4)
      iclocal(1) = (2.0d0)*(CON_KCALTOK)
      iclocal(3) = (180.0d0)*(onepi/180.0d0)
      iclocal(2) = 2.0d0
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     X-X-CA-HC [new compared to weiner1984]
      isame = isame + 1
      ninlocal(isame,1) = 'wild'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'CA'
      ninlocal(isame,4) = 'HC'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return

c     ******************************************************************
c     * special hydrogen bonding terms                                 *
c     ******************************************************************

      type = type + 1
      hbclocal(1) = (7557.0d0)*CON_KCALTOK
      hbclocal(2) = (2385.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H'
      hbnlocal(2) = 'NB'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (10238.0d0)*CON_KCALTOK
      hbclocal(2) = (3071.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H'
      hbnlocal(2) = 'NC'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (4019.0d0)*CON_KCALTOK
      hbclocal(2) = (1409.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H'
      hbnlocal(2) = 'O2'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (7557.0d0)*CON_KCALTOK
      hbclocal(2) = (2385.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H'
      hbnlocal(2) = 'O'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (7557.0d0)*CON_KCALTOK
      hbclocal(2) = (2385.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H'
      hbnlocal(2) = 'OH'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
c     --- slightly different from weiner1984
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H'
      hbnlocal(2) = 'S'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
c     --- slightly different from weiner1984
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H'
      hbnlocal(2) = 'SH'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (7557.0d0)*CON_KCALTOK
      hbclocal(2) = (2385.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HO'
      hbnlocal(2) = 'NB'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
c     --- new compared to weiner1984
      hbclocal(1) = (7557.0d0)*CON_KCALTOK
      hbclocal(2) = (2385.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HO'
      hbnlocal(2) = 'NC'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (4019.0d0)*CON_KCALTOK
      hbclocal(2) = (1409.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HO'
      hbnlocal(2) = 'O2'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (7557.0d0)*CON_KCALTOK
      hbclocal(2) = (2385.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HO'
      hbnlocal(2) = 'O'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (7557.0d0)*CON_KCALTOK
      hbclocal(2) = (2385.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HO'
      hbnlocal(2) = 'OH'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
c     --- slightly different from weiner1984
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HO'
      hbnlocal(2) = 'S'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
c     --- slightly different from weiner1984
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HO'
      hbnlocal(2) = 'SH'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (4019.0d0)*CON_KCALTOK
      hbclocal(2) = (1409.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H2'
      hbnlocal(2) = 'NB'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
c     --- new compared to weiner1984
      hbclocal(1) = (4019.0d0)*CON_KCALTOK
      hbclocal(2) = (1409.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H2'
      hbnlocal(2) = 'NC'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (4019.0d0)*CON_KCALTOK
      hbclocal(2) = (1409.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H2'
      hbnlocal(2) = 'O2'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (10238.0d0)*CON_KCALTOK
      hbclocal(2) = (3071.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H2'
      hbnlocal(2) = 'O'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (4019.0d0)*CON_KCALTOK
      hbclocal(2) = (1409.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H2'
      hbnlocal(2) = 'OH'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
c     --- slightly different from weiner1984
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H2'
      hbnlocal(2) = 'S'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
c     --- slightly different from weiner1984
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H2'
      hbnlocal(2) = 'SH'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (4019.0d0)*CON_KCALTOK
      hbclocal(2) = (1409.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H3'
      hbnlocal(2) = 'NB'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (4019.0d0)*CON_KCALTOK
      hbclocal(2) = (1409.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H3'
      hbnlocal(2) = 'NC'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (4019.0d0)*CON_KCALTOK
      hbclocal(2) = (1409.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H3'
      hbnlocal(2) = 'O2'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
c     --- very different values from weiner1984
      hbclocal(1) = (7557.0d0)*CON_KCALTOK
      hbclocal(2) = (2385.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H3'
      hbnlocal(2) = 'O'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
c     --- very different values from weiner1984
      hbclocal(1) = (7557.0d0)*CON_KCALTOK
      hbclocal(2) = (2385.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H3'
      hbnlocal(2) = 'OH'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
c     --- very different values from weiner1984
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H3'
      hbnlocal(2) = 'S'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
c     --- very different values from weiner1984
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'H3'
      hbnlocal(2) = 'SH'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (14184.0d0)*CON_KCALTOK
      hbclocal(2) = (3082.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HS'
      hbnlocal(2) = 'NB'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
c     --- new compared to weiner1984
      hbclocal(1) = (14184.0d0)*CON_KCALTOK
      hbclocal(2) = (3082.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HS'
      hbnlocal(2) = 'NC'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
c     --- very different values compared to weiner1984
      hbclocal(1) = (14184.0d0)*CON_KCALTOK
      hbclocal(2) = (3082.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HS'
      hbnlocal(2) = 'O2'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (14184.0d0)*CON_KCALTOK
      hbclocal(2) = (3082.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HS'
      hbnlocal(2) = 'O'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (14184.0d0)*CON_KCALTOK
      hbclocal(2) = (3082.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HS'
      hbnlocal(2) = 'OH'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
c     --- slightly different parameters from weiner1984
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HS'
      hbnlocal(2) = 'S'
      call twh_set_hbondnames(type,hbnlocal)

      type = type + 1
      hbclocal(1) = (265720.0d0)*CON_KCALTOK
c     --- slightly different parameters from weiner1984
      hbclocal(2) = (35029.0d0)*CON_KCALTOK
      hbclocal(3) = 0.5d0
      call twh_set_hbondcoeff(type,hbclocal)
      hbnlocal(1) = 'HS'
      hbnlocal(2) = 'SH'
      call twh_set_hbondnames(type,hbnlocal)

      nhbondtypes = type
      call twh_ffcheck(lfailure,8,type)
      if ( lfailure ) return

      return
      end
