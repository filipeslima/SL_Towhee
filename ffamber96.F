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
      subroutine twh_ffamber96(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the amber96      *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-20-2002 by M.G. Martin  *
c     * last modified 07-21-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
c     --- local variables
      integer ii,type
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal,ntslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
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
      double precision onepi

      ffname = 'Amber96'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

      type = 0
c     **************************************************************
c     * Amber param96 force field                *
c     **************************************************************
c     There are three version of what has become known as the
c     Cornell version of the Amber force field (as opposed to an 
c     earlier effort by Weiner et al. which is not currently coded 
c     up in Towhee).  The original param94 version (corn1995) is the
c     basis for the other versions.  param96 has a few torsions that
c     have different parameters.  The reference for that work is in 
c     a slightly obscure book, but the parameters are all listed on 
c     the Amber web site 
c     http://www.amber.ucsf.edu:80/amber/ff94.html
c     param99 was just released and has not yet been implemented 
c     into Towhee, although that work is planned for the near future
c     references:
c     [corn1995] W.D. Cornell, P. Cieplak, C.I. Bayly, I.R. Gould,
c     ---------- K.M. Merz Jr., D.M. Ferguson, D.C. Spellmeyer,
c     ---------- T. Fox, J.W. Caldwell, P.A. Kollman
c     ---------- `A second generation force field for the simulation
c     ---------- of proteins, nucleic acids, and organic molecules'
c     ---------- J. Am. Chem. Soc. 117, 5179-5197 (1995)
      
c     (C,CA,CM,CC,CV,CW,CR,CB,C*,CN,CK,CQ) any sp2 hybrid carbon
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CV'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CW'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CR'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CB'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C*'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CN'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CK'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0860d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CQ'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (CT) any sp3 hybrid carbon
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9080d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1094d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Cs) Cesium
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.3950d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0000806d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 132.90d0
      elelocal = 'Cs'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cs'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (F) Fluorine
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.75d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.061d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 18.998d0
      elelocal = ' F'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='F'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (H) H attached to N
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.600d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0157d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (H1) H attached to aliphatic C with 1 electron withdrawing
c     ---- group
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3870d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0157d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (H2) H attached to aliphatic C with 2 electron withdrawing
c     ---- groups
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.2870d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0157d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (H3) H attached to aliphatic C with 3 electron withdrawing
c     ---- groups
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.1870d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0157d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (H4) H attached to aromatic C with 1 electronegative 
c     ---- neighbor
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.4090d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0150d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (H5) H attached to aromatic C with 2 electronegative
c     ---- neighbors
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3590d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0150d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (HA) H attached aromatic C with no electronegative 
c     ---- neighbors
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.4590d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0150d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (HC) H attached to aliphatic C with no electron withdrawing
c     ---- groups
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.4870d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0157d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.008d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (HO) H attached to O in alcohols and acids
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.0d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.008d0
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      nbnlocal ='HO'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 


c     (HW) H attached to O in TIP3P water
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.0d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.008d0
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      nbnlocal ='HW'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (HP) H attached to C directly bonded to formally positive
c     ---- atom
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.100d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0157d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.008d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HP'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (HS) H attached to S
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.600d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0157d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.008d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (LP) lone pairs
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8680d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.00277d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 0.0d0
      elelocal = '--'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='LP'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (K) Potasium
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.6580d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.000328d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 39.098d0
      elelocal = ' K'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='K'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Li) Lithium
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.1370d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0183d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 6.941d0
      elelocal = 'Li'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Li'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (N,NA,NB,NC,N*,N2) any sp2 Nitrogen

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8240d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8240d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='NA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8240d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='NB'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8240d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='NC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8240d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N*'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8240d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (N3) sp3 Nitrogen
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.875d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (O,O2) sp2 Oxygen in amides or in anionic acids
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6612d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2100d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6612d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2100d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (OH) sp3 Oxygen in alcohols, tyrosine and protonated 
c     ---- carboxylic acids
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7210d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2104d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OH'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (OS) sp3 Oxygen in ethers
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6837d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (OW) sp3 Oxygen in TIP3P water
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7683d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1520d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OW'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (P) phosphorus in phosphates
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.100d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 30.974d0
      elelocal = ' P'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='P'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Rb) Rubidium
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.9560d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.00017d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 85.47d0
      elelocal = 'Rb'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Rb'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (S,SH) Sulfur in methionine and cysteine

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.00d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.25d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='S'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.00d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.25d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='SH'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Cl-) ionic chlorine
c     --- taken from the amber_param96.dat file
c     --- there they cite Smith and Dang JCP 1994,100:5,3757
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.47d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.10d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 35.45d0
      elelocal = 'Cl'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cl-'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return
      
c     ******************************************************************
c     * Amber param96 bond potentials                                  *
c     ******************************************************************
c     -- W.D. Cornell et al. J. Am. Chem. Soc. (1995), 117, 5179-5197

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.409d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - CA
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.419d0
      vclocal(1) = (447.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- C - CB
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CB'
c     --- CB - CN
      vnlocal(2,1) = 'CB'
      vnlocal(2,2) = 'CN'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.444d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - CM
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - CT
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (490.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - N
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'N'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.383d0
      vclocal(1) = (424.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - N*
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'N*'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.388d0
      vclocal(1) = (418.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - NA
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.358d0
      vclocal(1) = (457.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - NC
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'NC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.229d0
      vclocal(1) = (570.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - O
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'O'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.250d0
      vclocal(1) = (656.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - O2
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'O2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.364d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - OH
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'OH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.459d0
      vclocal(1) = (388.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C* - CB
      vnlocal(1,1) = 'C*'
      vnlocal(1,2) = 'CB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.495d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C* - CT
      vnlocal(1,1) = 'C*'
      vnlocal(1,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.352d0
      vclocal(1) = (546.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C* - CW
      vnlocal(1,1) = 'C*'
      vnlocal(1,2) = 'CW'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.080d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 11
c     --- CA - HA
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'H4'
c     --- CA - H4
      vnlocal(2,1) = 'CA'
      vnlocal(2,2) = 'HA'
c     --- C* - HC
      vnlocal(3,1) = 'C*'
      vnlocal(3,2) = 'HC'
c     --- CK - H5
      vnlocal(4,1) = 'CK'
      vnlocal(4,2) = 'H5'
c     --- CM - H4
      vnlocal(5,1) = 'CM'
      vnlocal(5,2) = 'H4'
c     --- CM - H5
      vnlocal(6,1) = 'CM'
      vnlocal(6,2) = 'H5'
c     --- CM - HA
      vnlocal(7,1) = 'CM'
      vnlocal(7,2) = 'HA'
c     --- CQ - H5
      vnlocal(8,1) = 'CQ'
      vnlocal(8,2) = 'H5'
c     --- CR - H5
      vnlocal(9,1) = 'CR'
      vnlocal(9,2) = 'H5'
c     --- CV - H4
      vnlocal(10,1) = 'CV'
      vnlocal(10,2) = 'H4'
c     --- CW - H4
      vnlocal(11,1) = 'CW'
      vnlocal(11,2) = 'H4'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.400d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CA - CA
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CA'
c     --- CA - CN
      vnlocal(2,1) = 'CA'
      vnlocal(2,2) = 'CN'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.404d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CA - CB
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.433d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CA - CM
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.510d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CA - CT
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CT'
c     --- CM - CT
      vnlocal(2,1) = 'CM'
      vnlocal(2,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.340d0
      vclocal(1) = (481.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CA - N2
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'N2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.381d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CW - NA
      vnlocal(1,1) = 'CW'
      vnlocal(1,2) = 'NA'
c     --- CA - NA
      vnlocal(2,1) = 'CA'
      vnlocal(2,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.339d0
      vclocal(1) = (483.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CA - NC
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'NC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.370d0
      vclocal(1) = (520.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CB - CB
      vnlocal(1,1) = 'CB'
      vnlocal(1,2) = 'CB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.374d0
      vclocal(1) = (436.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CB - N*
      vnlocal(1,1) = 'CB'
      vnlocal(1,2) = 'N*'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.391d0
      vclocal(1) = (414.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CB - NB
      vnlocal(1,1) = 'CB'
      vnlocal(1,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.354d0
      vclocal(1) = (461.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CB - NC
      vnlocal(1,1) = 'CB'
      vnlocal(1,2) = 'NC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.504d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CC - CT
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.375d0
      vclocal(1) = (512.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CC - CV
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'CV'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (518.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CC - CW
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'CW'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.385d0
      vclocal(1) = (422.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CC - NA
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.394d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CC - NB
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'NB'
c     --- CV - NB
      vnlocal(2,1) = 'CV'
      vnlocal(2,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (440.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CK - N*
      vnlocal(1,1) = 'CK'
      vnlocal(1,2) = 'N*'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.304d0
      vclocal(1) = (529.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CK - NB
      vnlocal(1,1) = 'CK'
      vnlocal(1,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.350d0
      vclocal(1) = (549.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CM - CM
      vnlocal(1,1) = 'CM'
      vnlocal(1,2) = 'CM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (448.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CM - N*
      vnlocal(1,1) = 'CM'
      vnlocal(1,2) = 'N*'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.380d0
      vclocal(1) = (428.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CN - NA
      vnlocal(1,1) = 'CN'
      vnlocal(1,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.324d0
      vclocal(1) = (502.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CQ - NC
      vnlocal(1,1) = 'CQ'
      vnlocal(1,2) = 'NC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.343d0
      vclocal(1) = (477.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CR - NA
      vnlocal(1,1) = 'CR'
      vnlocal(1,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (488.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CR - NB
      vnlocal(1,1) = 'CR'
      vnlocal(1,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (310.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - CT
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.380d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - F
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'F'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.090d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 5
c     --- CT - H1
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'H1'
c     --- CT - H2
      vnlocal(2,1) = 'CT'
      vnlocal(2,2) = 'H2'
c     --- CT - H3
      vnlocal(3,1) = 'CT'
      vnlocal(3,2) = 'H3'
c     --- CT - HC
      vnlocal(4,1) = 'CT'
      vnlocal(4,2) = 'HC'
c     --- CT - HP
      vnlocal(5,1) = 'CT'
      vnlocal(5,2) = 'HP'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.449d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CT - N
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'N'
c     --- CH - N
      vnlocal(2,1) = 'CT'
      vnlocal(2,2) = 'N'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.475d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - N*
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'N*'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.463d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - N2
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'N2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.471d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - N3
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'N3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.410d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CT - OH
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'OH'
c     --- CT - OS
      vnlocal(2,1) = 'CT'
      vnlocal(2,2) = 'OS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.810d0
      vclocal(1) = (227.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - S
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'S'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.810d0
      vclocal(1) = (237.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - SH
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'SH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.010d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 5
c     --- H - N
      vnlocal(1,1) = 'H'
      vnlocal(1,2) = 'N'
c     --- H - N*
      vnlocal(2,1) = 'H'
      vnlocal(2,2) = 'N*'
c     --- H - N2
      vnlocal(3,1) = 'H'
      vnlocal(3,2) = 'N2'
c     --- H - N3
      vnlocal(4,1) = 'H'
      vnlocal(4,2) = 'N3'
c     --- H - NA
      vnlocal(5,1) = 'H'
      vnlocal(5,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.960d0
      vclocal(1) = (553.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- HO - OH
      vnlocal(1,1) = 'HO'
      vnlocal(1,2) = 'OH'
c     --- HO - OS
      vnlocal(2,1) = 'HO'
      vnlocal(2,2) = 'OS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.336d0
      vclocal(1) = (274.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HS - SH
      vnlocal(1,1) = 'HS'
      vnlocal(1,2) = 'SH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.480d0
      vclocal(1) = (525.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- O2 - P
      vnlocal(1,1) = 'O2'
      vnlocal(1,2) = 'P'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.610d0
      vclocal(1) = (230.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- OH - P 
      vnlocal(1,1) = 'OH'
      vnlocal(1,2) = 'P'
c     --- OS - P 
      vnlocal(2,1) = 'OS'
      vnlocal(2,2) = 'P'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9572d0
      vclocal(1) = (553.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OW - HW
      vnlocal(1,1) = 'OW'
      vnlocal(1,2) = 'HW'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.038d0
      vclocal(1) = (166.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- S - S
      vnlocal(1,1) = 'S'
      vnlocal(1,2) = 'S'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Amber96 simple harmonic                                        *
c     ******************************************************************
c     -- W.D. Cornell et al. J. Am. Chem. Soc. (1995), 117, 5179-5197

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- C - CA - CA
c     --- CA - C,CA - CA,CB,CN
      nbslocal = 7
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CA'
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'CA'
      bnlocal(3,1) = 'CA'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'CB'
      bnlocal(4,1) = 'CA'
      bnlocal(4,2) = 'C'
      bnlocal(4,3) = 'CN'
      bnlocal(5,1) = 'CA'
      bnlocal(5,2) = 'CA'
      bnlocal(5,3) = 'CA'
      bnlocal(6,1) = 'CA'
      bnlocal(6,2) = 'CA'
      bnlocal(6,3) = 'CB'
      bnlocal(7,1) = 'CA'
      bnlocal(7,2) = 'CA'
      bnlocal(7,3) = 'CN'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- C - CA - HA
c     --- C* - CW - H4
c     --- CA - CA - H4,HA
c     --- CA - N2 - H
c     --- CB - CA - H4,HA
c     --- CC - CV,CW - H4 
c     --- CN - CA - HA
c     --- H - N,N2 - H
c     --- H4 - CV - NB : H4 - CW - NA
c     --- H5 - CR - NA,NB
      nbslocal = 16
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'HA'
      bnlocal(2,1) = 'C*'
      bnlocal(2,2) = 'CW'
      bnlocal(2,3) = 'H4'
      bnlocal(3,1) = 'CA'
      bnlocal(3,2) = 'CA'
      bnlocal(3,3) = 'H4'
      bnlocal(4,1) = 'CA'
      bnlocal(4,2) = 'CA'
      bnlocal(4,3) = 'HA'
      bnlocal(5,1) = 'CA'
      bnlocal(5,2) = 'N2'
      bnlocal(5,3) = 'H'
      bnlocal(6,1) = 'CB'
      bnlocal(6,2) = 'CA'
      bnlocal(6,3) = 'H4'
      bnlocal(7,1) = 'CB'
      bnlocal(7,2) = 'CA'
      bnlocal(7,3) = 'HA'
      bnlocal(8,1) = 'CC'
      bnlocal(8,2) = 'CV'
      bnlocal(8,3) = 'H4'
      bnlocal(9,1) = 'CC'
      bnlocal(9,2) = 'CW'
      bnlocal(9,3) = 'H4'
      bnlocal(10,1) = 'CN'
      bnlocal(10,2) = 'CA'
      bnlocal(10,3) = 'HA'
      bnlocal(11,1) = 'H'
      bnlocal(11,2) = 'N'
      bnlocal(11,3) = 'H'
      bnlocal(12,1) = 'H'
      bnlocal(12,2) = 'N2'
      bnlocal(12,3) = 'H'
      bnlocal(13,1) = 'H4'
      bnlocal(13,2) = 'CV'
      bnlocal(13,3) = 'NB'
      bnlocal(14,1) = 'H4'
      bnlocal(14,2) = 'CW'
      bnlocal(14,3) = 'NA'
      bnlocal(15,1) = 'H5'
      bnlocal(15,2) = 'CR'
      bnlocal(15,3) = 'NA'
      bnlocal(16,1) = 'H5'
      bnlocal(16,2) = 'CR'
      bnlocal(16,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.20d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- C - CB - CB
      nbslocal = 1
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CB - NB
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- C - CM - CM
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.70d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'CM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.70d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - CM - CT
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'CT'
c     --- CM - CM - CT
      bnlocal(2,1) = 'CM'
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.70d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- C - CM - H4,HA
c     --- CM - CM - H4,HA
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'H4'
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'HA'
      bnlocal(3,1) = 'CM'
      bnlocal(3,2) = 'CM'
      bnlocal(3,3) = 'H4'
      bnlocal(4,1) = 'CM'
      bnlocal(4,2) = 'CM'
      bnlocal(4,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.10d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CT - CT
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- C - CT - H1,HC,HP
c     --- C* - CT - HC
c     --- CA - CT - HC
c     --- CC - CT - HC
c     --- CM - CT - HC
c     --- CT - CT - H1,H2,HC,HP,N*,OH,OS
c     --- CT - N3 - H
c     --- H1 - CT - N,N*,N2,OH,OS,S,SH
c     --- H2 - CT - N*,OS
C     --- HP - CT - N3
c     --- N* - CT - OS
c     --- CT - N3 - CT
      nbslocal = 27
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'H1'
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'HC'
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'CT'
      bnlocal(3,3) = 'HP'
      bnlocal(4,1) = 'C*'
      bnlocal(4,2) = 'CT'
      bnlocal(4,3) = 'HC'
      bnlocal(5,1) = 'CA'
      bnlocal(5,2) = 'CT'
      bnlocal(5,3) = 'HC'
      bnlocal(6,1) = 'CC'
      bnlocal(6,2) = 'CT'
      bnlocal(6,3) = 'HC'
      bnlocal(7,1) = 'CM'
      bnlocal(7,2) = 'CT'
      bnlocal(7,3) = 'HC'
      bnlocal(8,1) = 'CT'
      bnlocal(8,2) = 'CT'
      bnlocal(8,3) = 'H1'
      bnlocal(9,1) = 'CT'
      bnlocal(9,2) = 'CT'
      bnlocal(9,3) = 'H2'
      bnlocal(10,1) = 'CT'
      bnlocal(10,2) = 'CT'
      bnlocal(10,3) = 'HC'
      bnlocal(11,1) = 'CT'
      bnlocal(11,2) = 'CT'
      bnlocal(11,3) = 'HP'
      bnlocal(12,1) = 'CT'
      bnlocal(12,2) = 'CT'
      bnlocal(12,3) = 'N*'
      bnlocal(13,1) = 'CT'
      bnlocal(13,2) = 'CT'
      bnlocal(13,3) = 'OH'
      bnlocal(14,1) = 'CT'
      bnlocal(14,2) = 'CT'
      bnlocal(14,3) = 'OS'
      bnlocal(15,1) = 'CT'
      bnlocal(15,2) = 'N3'
      bnlocal(15,3) = 'H'
      bnlocal(16,1) = 'H1'
      bnlocal(16,2) = 'CT'
      bnlocal(16,3) = 'N'
      bnlocal(17,1) = 'H1'
      bnlocal(17,2) = 'CT'
      bnlocal(17,3) = 'N*'
      bnlocal(18,1) = 'H1'
      bnlocal(18,2) = 'CT'
      bnlocal(18,3) = 'N2'
      bnlocal(19,1) = 'H1'
      bnlocal(19,2) = 'CT'
      bnlocal(19,3) = 'OH'
      bnlocal(20,1) = 'H1'
      bnlocal(20,2) = 'CT'
      bnlocal(20,3) = 'OS'
      bnlocal(21,1) = 'H1'
      bnlocal(21,2) = 'CT'
      bnlocal(21,3) = 'S'
      bnlocal(22,1) = 'H1'
      bnlocal(22,2) = 'CT'
      bnlocal(22,3) = 'SH'
      bnlocal(23,1) = 'H2'
      bnlocal(23,2) = 'CT'
      bnlocal(23,3) = 'N*'
      bnlocal(24,1) = 'H2'
      bnlocal(24,2) = 'CT'
      bnlocal(24,3) = 'OS'
      bnlocal(25,1) = 'HP'
      bnlocal(25,2) = 'CT'
      bnlocal(25,3) = 'N3'
      bnlocal(26,1) = 'N*'
      bnlocal(26,2) = 'CT'
      bnlocal(26,3) = 'OS'
      bnlocal(27,1) = 'CT'
      bnlocal(27,2) = 'N3'
      bnlocal(27,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.10d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CT - N
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.20d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CT - N3
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'N3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.90d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - N - CT
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- C - N - H
c     --- CC - NA - H
c     --- CR - NA - H
c     --- CW - NA - H
      nbslocal = 4
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'H'
      bnlocal(2,1) = 'CC'
      bnlocal(2,2) = 'NA'
      bnlocal(2,3) = 'H'
      bnlocal(3,1) = 'CR'
      bnlocal(3,2) = 'NA'
      bnlocal(3,3) = 'H'
      bnlocal(4,1) = 'CW'
      bnlocal(4,2) = 'NA'
      bnlocal(4,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - N* - CM
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'CM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - N* - CT
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.20d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- C - N* - H
      nbslocal = 1
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - NA - C
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - NA - CA
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.80d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- C - NA - H
      nbslocal = 1
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - NC - CA
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NC'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.00d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - OH - HO
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'OH'
      bnlocal(1,3) = 'HO'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 134.90d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- C* - CB - CA
      nbslocal = 1
      bnlocal(1,1) = 'C*'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.80d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- C* - CB - CN
      nbslocal = 1
      bnlocal(1,1) = 'C*'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CN'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.60d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C* - CT - CT
      bnlocal(1,1) = 'C*'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.70d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C* - CW - NA
      bnlocal(1,1) = 'C*'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- CA - C - OH
c     --- CA - CA - CT
c     --- CC - CV - NB
c     --- CC- CW - NA 
c     --- CC - NA - CR
c     --- CR - NA - CW
c     --- CT - CC - CV,CW,NA,NB
c     --- CV - CC - NA 
c     --- CW - CC - NA,NB
c     --- N2 - CA - N2
c     --- NA - CR - NA,NB
      nbslocal = 16
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OH'
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'CT'
      bnlocal(3,1) = 'CC'
      bnlocal(3,2) = 'CV'
      bnlocal(3,3) = 'NB'
      bnlocal(4,1) = 'CC'
      bnlocal(4,2) = 'CW'
      bnlocal(4,3) = 'NA'
      bnlocal(5,1) = 'CC'
      bnlocal(5,2) = 'NA'
      bnlocal(5,3) = 'CR'
      bnlocal(6,1) = 'CR'
      bnlocal(6,2) = 'NA'
      bnlocal(6,3) = 'CW'
      bnlocal(7,1) = 'CT'
      bnlocal(7,2) = 'CC'
      bnlocal(7,3) = 'CV'
      bnlocal(8,1) = 'CT'
      bnlocal(8,2) = 'CC'
      bnlocal(8,3) = 'CW'
      bnlocal(9,1) = 'CT'
      bnlocal(9,2) = 'CC'
      bnlocal(9,3) = 'NA'
      bnlocal(10,1) = 'CT'
      bnlocal(10,2) = 'CC'
      bnlocal(10,3) = 'NB'
      bnlocal(11,1) = 'CV'
      bnlocal(11,2) = 'CC'
      bnlocal(11,3) = 'NA'
      bnlocal(12,1) = 'CW'
      bnlocal(12,2) = 'CC'
      bnlocal(12,3) = 'NA'
      bnlocal(13,1) = 'CW'
      bnlocal(13,2) = 'CC'
      bnlocal(13,3) = 'NB'
      bnlocal(14,1) = 'N2'
      bnlocal(14,2) = 'CA'
      bnlocal(14,3) = 'N2'
      bnlocal(15,1) = 'NA'
      bnlocal(15,2) = 'CR'
      bnlocal(15,3) = 'NA'
      bnlocal(16,1) = 'NA'
      bnlocal(16,2) = 'CR'
      bnlocal(16,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.30d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- CA - CB - CB 
      nbslocal = 1
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.20d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- CA - CB - CN 
      nbslocal = 1
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CN'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 132.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - CB - NB
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- CA - CM - CM
      nbslocal = 1
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'CM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.30d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CA - CM - H4,HA
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'H4'
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.70d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- CA - CN - CB
      nbslocal = 1
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CN'
      bnlocal(1,3) = 'CB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 132.80d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - CN - NA
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CN'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.00d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - CT - CT 
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.20d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - N2 - CT
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'N2'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- CA - NA - H
      nbslocal = 1
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - NC - CB
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'NC'
      bnlocal(1,3) = 'CB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.6d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CA - NC - CQ
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'NC'
      bnlocal(1,3) = 'CQ'
c     --- N* - C - NC
      bnlocal(2,1) = 'N*'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - C - NA
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.80d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - C - O
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - C* - CT
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'C*'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.40d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - C* - CW
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'C*'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.50d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CA - N2
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'N2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CA - NC
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CB - N*
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'N*'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CB - NB
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.70d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CB - NC
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CN - NA
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CN'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - N* - CK
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'CK'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.80d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - N* - CT
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.80d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- CB - N* - H
      nbslocal = 1
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.80d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - NB - CK
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'NB'
      bnlocal(1,3) = 'CK'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - NC - CQ
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'NC'
      bnlocal(1,3) = 'CQ'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CC - CT - CT
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.10d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CC'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CC - NB - CR
      bnlocal(1,1) = 'CC'
      bnlocal(1,2) = 'NB'
      bnlocal(1,3) = 'CR'
c     --- CR - NB - CV
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'NB'
      bnlocal(2,3) = 'CV'
c     --- CT - C - O2
      bnlocal(3,1) = 'CT'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'O2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.80d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CK - N* - CT
      bnlocal(1,1) = 'CK'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CK - N* - H
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.80d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CK'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.10d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CM - C - NA
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.30d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CM - C - O
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.10d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CM - CA - N2
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'N2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.50d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CM - CA - NC
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CM - CM - N* 
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'N*'
c     --- CM - N* - CT
      bnlocal(2,1) = 'CM'
      bnlocal(2,2) = 'N*'
      bnlocal(2,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CM - N* - H
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.20d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CN - NA - CW
      bnlocal(1,1) = 'CN'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CN - NA - H
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.10d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CN'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - C - N
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.40d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - C - O
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - C* - CW
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'C*'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (40.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - CT - CT
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.70d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - CT - N
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CT - CT - N2,N3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.20d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'N2'
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'N3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.70d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - CT - S
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'S'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.60d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - CT - SH
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'SH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - N - CT
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CT - N - H
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.04d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CT - N2 - H
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.40d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'N2'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.50d0
      bclocal(1) = (55.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - OH - HO
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'OH'
      bnlocal(1,3) = 'HO'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (60.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - OS - CT
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - OS - P
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'P'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 98.90d0
      bclocal(1) = (62.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - S - CT
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.70d0
      bclocal(1) = (68.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - S - S
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'S'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 96.00d0
      bclocal(1) = (43.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - SH - HS
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'SH'
      bnlocal(1,3) = 'HS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.10d0
      bclocal(1) = (77.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- F - CT - F
      bnlocal(1,1) = 'F'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'F'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 6
c     --- F - CT - H1
      bnlocal(1,1) = 'F'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'H1'
c     --- H - N3 - H 
      bnlocal(2,1) = 'H'
      bnlocal(2,2) = 'N3'
      bnlocal(2,3) = 'H'
c     --- H1 - CT - H1 
      bnlocal(3,1) = 'H1'
      bnlocal(3,2) = 'CT'
      bnlocal(3,3) = 'H1'
c     --- H2 - CT - H2
      bnlocal(4,1) = 'H2'
      bnlocal(4,2) = 'CT'
      bnlocal(4,3) = 'H2'
C     --- HC - CT - HC
      bnlocal(5,1) = 'HC'
      bnlocal(5,2) = 'CT'
      bnlocal(5,3) = 'HC'
C     --- HP - CT - HP
      bnlocal(6,1) = 'HP'
      bnlocal(6,2) = 'CT'
      bnlocal(6,3) = 'HP'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.10d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- H4 - CM - N*
      bnlocal(1,1) = 'H4'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'N*'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.05d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- H5 - CK - N*,NB
      bnlocal(1,1) = 'H5'
      bnlocal(1,2) = 'CK'
      bnlocal(1,3) = 'N*'
      bnlocal(2,1) = 'H5'
      bnlocal(2,2) = 'CK'
      bnlocal(2,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.45d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- H5 - CQ - NC
      bnlocal(1,1) = 'H5'
      bnlocal(1,2) = 'CQ'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.50d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- HO - OH - P
      bnlocal(1,1) = 'HO'
      bnlocal(1,2) = 'OH'
      bnlocal(1,3) = 'P'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 92.07d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- HS - SH - HS
      bnlocal(1,1) = 'HS'
      bnlocal(1,2) = 'SH'
      bnlocal(1,3) = 'HS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- HW - OW - HW
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.52d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HW'
      bnlocal(1,2) = 'OW'
      bnlocal(1,3) = 'HW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.90d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N - C - O
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- N* - C - NA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'N*'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- N* - C - O
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.90d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'N*'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N* - CB - NC
      bnlocal(1,1) = 'N*'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.90d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N* - CK - NB
      bnlocal(1,1) = 'N*'
      bnlocal(1,2) = 'CK'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N2 - CA - NA
      bnlocal(1,1) = 'N2'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N2 - CA - NC
      bnlocal(1,1) = 'N2'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.60d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NA - C - O
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NA - CA - NC
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.50d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NC - C - O
      bnlocal(1,1) = 'NC'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.10d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NC - CQ - NC
      bnlocal(1,1) = 'NC'
      bnlocal(1,2) = 'CQ'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.00d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- O - C - O
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
c     --- O2 - C - O2
      bnlocal(2,1) = 'O2'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'O2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.90d0
      bclocal(1) = (140.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- O2 - P - O2
      bnlocal(1,1) = 'O2'
      bnlocal(1,2) = 'P'
      bnlocal(1,3) = 'O2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.23d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- O2 - P - OH
      bnlocal(1,1) = 'O2'
      bnlocal(1,2) = 'P'
      bnlocal(1,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.23d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- O2 - P - OS
      bnlocal(1,1) = 'O2'
      bnlocal(1,2) = 'P'
      bnlocal(1,3) = 'OS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 102.60d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- OH - P - OS
      bnlocal(1,1) = 'OH'
      bnlocal(1,2) = 'P'
      bnlocal(1,3) = 'OS'
c     --- OS - P - OS
      bnlocal(2,1) = 'OS'
      bnlocal(2,2) = 'P'
      bnlocal(2,3) = 'OS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- P - OS - P
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'P'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'P'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Amber96 standard dihedral                                      *
c     ******************************************************************
c     --- These torsional parameters are taken from the param96.dat
c     --- file listed on the Amber web site.  They are slightly 
c     --- different than the parameters in the Cornell et al. paper

c     --- get pi
      call twh_constant_pi(GLB_GET,onepi)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((14.50d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- X - C - CA - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CA'
      tnlocal(1,4) = 'wild'
c     --- X - CA - CA - X
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CA'
      tnlocal(2,3) = 'CA'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((12.00d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C - CB - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((8.70d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C - CM - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CM'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((5.80d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C - N* - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'N*'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((5.40d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C - NA - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'NA'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((8.00d0)/(2.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C - NC - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'NC'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((1.80d0)/(2.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C - OH - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'OH'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

C     --- Note: this torsion has no energy contribution from the cosine
c     --- term it is still important to the logic of the code (and the 
c     --- 1-4 term) to include this torsion in the input file.
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 8
c     --- X - C - CT - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'wild'
c     --- X - CA - CT - X
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CA'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'wild'
c     --- X - CM - CT - X
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CM'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'wild'
c     --- X - CT - N - X
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'N'
      tnlocal(4,4) = 'wild'
c     --- X - CT - N* - X
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'CT'
      tnlocal(5,3) = 'N*'
      tnlocal(5,4) = 'wild'
c     --- X - CT - N2 - X
      tnlocal(6,1) = 'wild'
      tnlocal(6,2) = 'CT'
      tnlocal(6,3) = 'N2'
      tnlocal(6,4) = 'wild'
c     --- X - C* - CT - X
      tnlocal(7,1) = 'wild'
      tnlocal(7,2) = 'C*'
      tnlocal(7,3) = 'CT'
      tnlocal(7,4) = 'wild'
c     --- X - CC - CT - X
      tnlocal(8,1) = 'wild'
      tnlocal(8,2) = 'CC'
      tnlocal(8,3) = 'CT'
      tnlocal(8,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((14.00d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CA - CB - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'CB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((10.20d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CA - CM - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'CM'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((9.60d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CA - N2 - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'N2'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((6.00d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CA - NA - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'NA'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((9.60d0)/(2.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CA - NC - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'NC'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((21.80d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CB - CB - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CB'
      tnlocal(1,3) = 'CB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((6.60d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CB - N* - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CB'
      tnlocal(1,3) = 'N*'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((5.10d0)/(2.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CB - NB - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CB'
      tnlocal(1,3) = 'NB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((8.30d0)/(2.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CB - NC - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CB'
      tnlocal(1,3) = 'NC'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((6.80d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CK - N* - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CK'
      tnlocal(1,3) = 'N*'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((20.00d0)/(2.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CK - NB - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CK'
      tnlocal(1,3) = 'NB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((26.60d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CM - CM - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CM'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((7.40d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CM - N* - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'N*'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((13.60d0)/(2.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CQ - NC - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CQ'
      tnlocal(1,3) = 'NC'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- be careful as there are some torsions below which also have 
c     --- ct-ct as the center interaction (with oxygen groups on the 
c     --- ends)  Completely defined torsions take precedence over wild
c     --- cards.
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((1.40d0)/(9.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- X - CT - CT - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'wild'
c     --- X - CT - N3 - X
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'N3'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((0.50d0)/(3.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CT - OH - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OH'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((1.15d0)/(3.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CT - OS - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((0.75d0)/(3.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- X - OH - P - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'OH'
      tnlocal(1,3) = 'P'
      tnlocal(1,4) = 'wild'
c     --- X - OS - P - X
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'OS'
      tnlocal(2,3) = 'P'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((10.0d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C - N - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((1.00d0)/(3.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CT - S - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((0.75d0)/(3.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CT - SH - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'SH'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((6.70d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C* - CB - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C*'
      tnlocal(1,3) = 'CB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((26.10d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C* - CW - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C*'
      tnlocal(1,3) = 'CW'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((14.50d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CA - CN - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'CN'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((12.00d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CB - CN - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CB'
      tnlocal(1,3) = 'CN'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((20.60d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CC - CV - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CC'
      tnlocal(1,3) = 'CV'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((21.50d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CC - CW - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CC'
      tnlocal(1,3) = 'CW'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((5.60d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CC - NA - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CC'
      tnlocal(1,3) = 'NA'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((4.80d0)/(2.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- X - CC - NB - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CC'
      tnlocal(1,3) = 'NB'
      tnlocal(1,4) = 'wild'
c     --- X - CV - NB - X
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CV'
      tnlocal(2,3) = 'NB'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((6.10d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CN - NA - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CN'
      tnlocal(1,3) = 'NA'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((9.30d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CR - NA - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CR'
      tnlocal(1,3) = 'NA'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((10.00d0)/(2.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CR - NB - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CR'
      tnlocal(1,3) = 'NB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = ((6.00d0)/(4.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CW - NA - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CW'
      tnlocal(1,3) = 'NA'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = ((0.383d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(4) = ((0.1d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      tclocal(5) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - OS - CT
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = ((0.30d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(4) = ((0.85d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(5) = 1.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- C - N - CT - C
      tnlocal(1,1) = 'C'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C'
c     --- N - CT - C - N
      tnlocal(2,1) = 'N'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'C'
      tnlocal(2,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = ((0.50d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 4.0d0
      tclocal(4) = ((0.15d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      tclocal(5) = 3.0d0
      tclocal(7) = ((0.53d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      tclocal(8) = 1.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - N - C
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'C'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = ((0.100d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 4.0d0
      tclocal(4) = ((0.07d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(5) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - C - N
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = ((2.50d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(4) = ((2.00d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(5) = 1.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- H - N - C - O
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = ((3.50d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(4) = ((0.60d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(5) = 3.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - S - S - CT
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'S'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = ((0.144d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(4) = ((1.00d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(5) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- OS - CT - CT- OS
      tnlocal(1,1) = 'OS'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OS'
c     --- OS - CT - CT - OH
      tnlocal(2,1) = 'OS'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'OH'
c     --- OH - CT - CT - OH
      tnlocal(3,1) = 'OH'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'OH'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = ((0.25d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(2) = 3.0d0
      tclocal(4) = ((1.20d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(5) = 2.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- OH - P - OS - CT
      tnlocal(1,1) = 'OH'
      tnlocal(1,2) = 'P'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT'
c     --- OS - P - OS - CT
      tnlocal(2,1) = 'OS'
      tnlocal(2,2) = 'P'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = ((0.50d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(4) = ((2.50d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(5) = 1.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- OS - CT - N* - CK
      tnlocal(1,1) = 'OS'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N*'
      tnlocal(1,4) = 'CK'
c     --- OS - CT - N* - CM
      tnlocal(2,1) = 'OS'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'N*'
      tnlocal(2,4) = 'CM'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 5.0d0/6.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = ((0.50d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(2) = 2.0d0
      tclocal(4) = ((2.50d0)/(1.0d0))*(CON_KCALTOK)
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(5) = 1.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- OS - CT - N* - CK
      tnlocal(1,1) = 'OS'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N*'
      tnlocal(1,4) = 'CK'
c     --- OS - CT - N* - CM
      tnlocal(2,1) = 'OS'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'N*'
      tnlocal(2,4) = 'CM'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Amber96 improper dihedral                                      *
c     * Amber improper parameter decoder                               *
c     * impcoeff(1)*[ 1 + Cos( impcoeff(2)* theta - impcoeff(3)mm) ]   *
c     ******************************************************************
C     --- These torsional parameters are taken from the param96.dat
c     --- file listed on the Amber web site.  They are slightly 
c     --- different than the parameters in the Cornell et al. paper
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
      nislocal = 4
c     --- X - X - C - O
      ninlocal(1,1) = 'wild'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'C'
      ninlocal(1,4) = 'O'
c     --- X - O2 - C - O2
      ninlocal(2,1) = 'wild'
      ninlocal(2,2) = 'O2'
      ninlocal(2,3) = 'C'
      ninlocal(2,4) = 'O2'
c     --- X - N2 - CA - N2
      ninlocal(3,1) = 'wild'
      ninlocal(3,2) = 'N2'
      ninlocal(3,3) = 'CA'
      ninlocal(3,4) = 'N2'
c     --- CT - O - C - OH
      ninlocal(4,1) = 'CT'
      ninlocal(4,2) = 'O'
      ninlocal(4,3) = 'C'
      ninlocal(4,4) = 'OH'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,4)
      iclocal(1) = (1.0d0)*(CON_KCALTOK)
      iclocal(3) = (180.0d0)*(onepi/180.0d0)
      iclocal(2) = 2.0d0
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 6
c     --- X - X - N - H
      ninlocal(1,1) = 'wild'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'N'
      ninlocal(1,4) = 'H'
c     --- X - X - N2 - H
      ninlocal(2,1) = 'wild'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'N2'
      ninlocal(2,4) = 'H'
c     --- X - X - NA - H
      ninlocal(3,1) = 'wild'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'NA'
      ninlocal(3,4) = 'H'
c     --- X - CT - N - CT
      ninlocal(4,1) = 'wild'
      ninlocal(4,2) = 'CT'
      ninlocal(4,3) = 'N'
      ninlocal(4,4) = 'CT'
c     --- CK - CB - N* - CT
      ninlocal(5,1) = 'CK'
      ninlocal(5,2) = 'CB'
      ninlocal(5,3) = 'N*'
      ninlocal(5,4) = 'CT'
c     --- CM - C - N* - CT
      ninlocal(6,1) = 'CM'
      ninlocal(6,2) = 'C'
      ninlocal(6,3) = 'N*'
      ninlocal(6,4) = 'CT'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,4)
      iclocal(1) = (1.1d0)*(CON_KCALTOK)
      iclocal(3) = (180.0d0)*(onepi/180.0d0)
      iclocal(2) = 2.0d0
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 21
c     --- X - X - CA - HA
      ninlocal(1,1) = 'wild'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'CA'
      ninlocal(1,4) = 'HA'
c     --- X - X - CW - H4
      ninlocal(2,1) = 'wild'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'CW'
      ninlocal(2,4) = 'H4'
c     --- X - X - CR - H5
      ninlocal(3,1) = 'wild'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'CR'
      ninlocal(3,4) = 'H5'
c     --- X - X - CV - H4
      ninlocal(4,1) = 'wild'
      ninlocal(4,2) = 'wild'
      ninlocal(4,3) = 'CV'
      ninlocal(4,4) = 'H4'
c     --- X - X - CQ - H5
      ninlocal(5,1) = 'wild'
      ninlocal(5,2) = 'wild'
      ninlocal(5,3) = 'CQ'
      ninlocal(5,4) = 'H5'
c     --- X - X - CK - H5
      ninlocal(6,1) = 'wild'
      ninlocal(6,2) = 'wild'
      ninlocal(6,3) = 'CK'
      ninlocal(6,4) = 'H5'
c     --- X - X - CM - H4
      ninlocal(7,1) = 'wild'
      ninlocal(7,2) = 'wild'
      ninlocal(7,3) = 'CM'
      ninlocal(7,4) = 'H4'
c     --- X - X - CM - HA
      ninlocal(8,1) = 'wild'
      ninlocal(8,2) = 'wild'
      ninlocal(8,3) = 'CM'
      ninlocal(8,4) = 'HA'
c     --- X - X - CA - H4
      ninlocal(9,1) = 'wild'
      ninlocal(9,2) = 'wild'
      ninlocal(9,3) = 'CA'
      ninlocal(9,4) = 'H4'
c     --- X - X - CA - H5
      ninlocal(10,1) = 'wild'
      ninlocal(10,2) = 'wild'
      ninlocal(10,3) = 'CA'
      ninlocal(10,4) = 'H5'
c     --- CM - C - CM - CT
      ninlocal(11,1) = 'CM'
      ninlocal(11,2) = 'C'
      ninlocal(11,3) = 'CM'
      ninlocal(11,4) = 'CT'
c     --- NA - CV - CC - CT
      ninlocal(12,1) = 'NA'
      ninlocal(12,2) = 'CV'
      ninlocal(12,3) = 'CC'
      ninlocal(12,4) = 'CT'
c     --- NA - CW - CC - CT
      ninlocal(13,1) = 'NA'
      ninlocal(13,2) = 'CW'
      ninlocal(13,3) = 'CC'
      ninlocal(13,4) = 'CT'
c     --- NB - CW - CC - CT
      ninlocal(14,1) = 'NB'
      ninlocal(14,2) = 'CW'
      ninlocal(14,3) = 'CC'
      ninlocal(14,4) = 'CT'
c     --- CW - CB - C* - CT
      ninlocal(15,1) = 'CW'
      ninlocal(15,2) = 'CB'
      ninlocal(15,3) = 'C*'
      ninlocal(15,4) = 'CT'
c     --- CA - CA - CA - CT
      ninlocal(16,1) = 'CA'
      ninlocal(16,2) = 'CA'
      ninlocal(16,3) = 'CA'
      ninlocal(16,4) = 'CT'
c     --- C - CM - CM - CT
      ninlocal(17,1) = 'C'
      ninlocal(17,2) = 'CM'
      ninlocal(17,3) = 'CM'
      ninlocal(17,4) = 'CT'
c     --- NC - CM - CA - N2
      ninlocal(18,1) = 'NC'
      ninlocal(18,2) = 'CM'
      ninlocal(18,3) = 'CA'
      ninlocal(18,4) = 'N2'
c     --- CB - NC - CA - N2
      ninlocal(19,1) = 'CB'
      ninlocal(19,2) = 'NC'
      ninlocal(19,3) = 'CA'
      ninlocal(19,4) = 'N2'
c     --- NA - NC - CA - N2
      ninlocal(20,1) = 'NA'
      ninlocal(20,2) = 'NC'
      ninlocal(20,3) = 'CA'
      ninlocal(20,4) = 'N2'
c     --- CA - CA - C - OH
      ninlocal(21,1) = 'CA'
      ninlocal(21,2) = 'CA'
      ninlocal(21,3) = 'C'
      ninlocal(21,4) = 'OH'
      call twh_set_impnames(type,nislocal,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return

      return
      end
