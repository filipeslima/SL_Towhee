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
      subroutine twh_ffcharmm22(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the charmm22     *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-20-2002 by M.G. Martin  *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
c     --- local variables
      integer ii,type,isame
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

      ffname = 'Charmm22'
      type = 0
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      
c     **************************************************************
c     * Charmm22 force field                                       *
c     **************************************************************
c     NOTE: there are more charmm22 parameters listed down below
c     starting at type 141.  These were added after the original
c     set as they are used for the schiff base in rhodopsin
c     (M.G. Martin 08-21-2001)
c     references:
c     [mack1998] A.D. MacKerell Jr., D. Bashford, M. Bellott, 
c     ---------- R.L. Dunbrack Jr., J.D. Evanseck, M.J. Field, 
c     ---------- S. Fischer, J. Gao, H. Guo, S. Ha, 
c     ---------- D. Joseph-McCarthy, L. Kuchnir, K. Kuczera, 
c     ---------- F.T.K. Lau, C. Mattos, S. Michnick, T. Ngo, 
c     ---------- D.T. Nguyen, B. Prodhom, W.E. Reiher III, B. Roux,
c     ---------- M. Schlenkrich, J.C. Smith, R. Stote, J. Straub, 
c     ---------- M. Watanabe, J. Wiorkiewicz-Kuczera, D. Yin, 
c     ---------- M. Karplus
c     ---------- `All-atom empirical potential for molecular 
c     ---------- modeling and dynamics studies of proteins'
c     ---------- J. Phys. Chem. B 102, 3586-3616 (1998) 
c     ---------- {see the supplemental info}
      
c     (C) polar C
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.0d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1100d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (CA) aromatic C
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9924d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
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
      nbclocal(1) = (2.00d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CC) carbonyl C for sidechains asn,asp,gln,glu
      nbnlocal ='CC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.00d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CD) carbonyl C for none amides, asp,glu,cter
      nbnlocal ='CD'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.10d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1100d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CM) heme CO carbon
      nbnlocal ='CM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.1750d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0550d0)*(CON_KCALTOK)
      nbclocal(3) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.0100d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CP2) tetrahedral C (proline CB/CG)
      nbnlocal ='CP2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.1750d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0550d0)*(CON_KCALTOK)
      nbclocal(3) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.0100d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CP3) tetrahedral C (proline CD)
      nbnlocal ='CP3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.1750d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0550d0)*(CON_KCALTOK)
      nbclocal(3) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.0100d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT2) sp3 hybrid carbon bonded to 2 hydrogens
      nbnlocal ='CT2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0900d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPA) heme alpha-C
      nbnlocal ='CPA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0900d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPB) heme beta-C
      nbnlocal ='CPB'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0900d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPM) heme meso-C
      nbnlocal ='CPM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0500d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPH1) his CG and CD2 carbons
      nbnlocal ='CPH1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0500d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPH2) his CE1 carbon
      nbnlocal ='CPH2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CPT) trp C between rings
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0900d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.0900d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CPT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CS) thiolate carbon
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.200d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1100d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.2750d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0200d0)*(CON_KCALTOK)
      nbclocal(3) = (1.90d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.01d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT1) sp3 hybrid carbon bonded to 1 hydrogens
      nbnlocal ='CT1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.2750d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0200d0)*(CON_KCALTOK)
      nbclocal(3) = (1.90d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.01d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CP1) tetrahedral C (proline CA)
      nbnlocal ='CP1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.0600d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0800d0)*(CON_KCALTOK)
      nbclocal(3) = (1.90d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.01d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT3) sp3 hybrid carbon bonded to 3 hydrogens
      nbnlocal ='CT3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (CY) trp C in pyrrole ring
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9924d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CY'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.2245d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0460d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (H) sp3 hybrid hydrogen bonded to polar group (polar H)
      nbnlocal ='H'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.2245d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0460d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HC) N-TER H
      nbnlocal ='HC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.2245d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0460d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HT) TIP3P water hydrogen
      nbnlocal ='HT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.32d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0220d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HA) hydrogen on an aliphatic chain (nonpolar H)
      nbnlocal ='HA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.32d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0220d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HB) backbone H
      nbnlocal ='HB'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (HP) aromatic H
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3582d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0300d0)*(CON_KCALTOK)
      nbclocal(3) = (1.3582d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.0300d0)*(CON_KCALTOK)
      masslocal = 1.00800d0
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

c         (HR1) his HE1, (+)his HG, HD2
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.9000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0460d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HR1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.7000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0460d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HR2) (+)his HE1
      nbnlocal ='HR2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c         (HR3) neutral his HG,HD2
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.4680d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0078d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HR3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c         (HS) thiol hydrogen
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.4500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
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
         
c         (N) proline N
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.0001d0)*(CON_KCALTOK)
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
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NC2) guanidinium nitrogen
      nbnlocal ='NC2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NH2) amide nitrogen
      nbnlocal ='NH2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NP) proline ring NH2+ (N-terminal)
      nbnlocal ='NP'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NPH) heme pyrrole N
      nbnlocal ='NPH'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NR1) neutral his protonated ring nitrogen
      nbnlocal ='NR1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NR2) neutral his unprotonated ring nitrogen
      nbnlocal ='NR2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NY) trp N in pyrrole ring
      nbnlocal ='NY'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NH3) ammonium nitrogen
      nbnlocal ='NH3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NR3) charged his ring nitrogen
      nbnlocal ='NR3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (NH1) peptide nitrogen
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.8500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = (1.5500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.2000d0)*(CON_KCALTOK)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='NH1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1200d0)*(CON_KCALTOK)
      nbclocal(3) = (1.400d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.1200d0)*(CON_KCALTOK)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (O) carbonyl oxygen
      nbnlocal ='O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1200d0)*(CON_KCALTOK)
      nbclocal(3) = (1.400d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.1200d0)*(CON_KCALTOK)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OB) carbonyl oxygen in acetic acid
      nbnlocal ='OB'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1200d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OC) carboxylate oxygen
      nbnlocal ='OC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1200d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OM) heme CO/O2 oxygen
      nbnlocal ='OM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.770d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.15210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OH1) hydroxyl oxygen
      nbnlocal ='OH1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.770d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.15210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OS) ester oxygen
      nbnlocal ='OS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7682d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1521d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OT) TIP3P water oxygen
      nbnlocal ='OT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.6500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 55.847d0
      elelocal = 'Fe'
      bplocal = 'heme'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (FE) heme iron
      nbnlocal ='FE'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (S) sulfur
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.4500d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
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

c     (SM) sulfur C-S-S-C type
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9750d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.3800d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='SM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (SS) thiolate sulfur 
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.2000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.4700d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='SS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Charmm22 bond potentials                                       *
c     ******************************************************************

c     --- C - C
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3350d0
      vclocal(1) = (600.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CA - CA
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3750d0
      vclocal(1) = (305.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CP1 - C,CC
c     --- CT1,CT2,CT3 - C
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4900d0
      vclocal(1) = (250.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 5
      vnlocal(1,1) = 'CP1'
      vnlocal(1,2) = 'C'
      vnlocal(2,1) = 'CP1'
      vnlocal(2,2) = 'CC'
      vnlocal(3,1) = 'CT1'
      vnlocal(3,2) = 'C'
      vnlocal(4,1) = 'CT2'
      vnlocal(4,2) = 'C'
      vnlocal(5,1) = 'CT3'
      vnlocal(5,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CP1 - CD
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4900d0
      vclocal(1) = (200.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CP1'
      vnlocal(1,2) = 'CD'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CP2 - CP1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5270d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CP2'
      vnlocal(1,2) = 'CP1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5370d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CP2 - CP2 [mack1998]
      vnlocal(1,1) = 'CP2'
      vnlocal(1,2) = 'CP2'
c     --- CP3 - CP2 [mack1998]
      vnlocal(2,1) = 'CP3'
      vnlocal(2,2) = 'CP2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3800d0
      vclocal(1) = (450.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CPB - C
      vnlocal(1,1) = 'CPB'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CPB - CPA
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4432d0
      vclocal(1) = (299.80d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CPB'
      vnlocal(1,2) = 'CPA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CPB - CPB
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3464d0
      vclocal(1) = (340.70d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CPB'
      vnlocal(1,2) = 'CPB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CPH1 - CPH1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3600d0
      vclocal(1) = (410.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CPH1'
      vnlocal(1,2) = 'CPH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CPM - CPA
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3716d0
      vclocal(1) = (360.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CPM'
      vnlocal(1,2) = 'CPA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CPT - CA
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3680d0
      vclocal(1) = (305.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CPT'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CPT - CPT
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.400d0
      vclocal(1) = (360.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CPT'
      vnlocal(1,2) = 'CPT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (200.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 6
c     --- CT1 - CC [mack1998]
      vnlocal(1,1) = 'CT1'
      vnlocal(1,2) = 'CC'
c     --- CT2 - CC [mack1998]
      vnlocal(2,1) = 'CT2'
      vnlocal(2,2) = 'CC'
c     --- CT3 - CC [mack1998]
      vnlocal(3,1) = 'CT3'
      vnlocal(3,2) = 'CC'
c     --- CT1 - CD [mack1998]
      vnlocal(4,1) = 'CT1'
      vnlocal(4,2) = 'CD'
c     --- CT2 - CD [mack1998]
      vnlocal(5,1) = 'CT2'
      vnlocal(5,2) = 'CD'
c     --- CT3 - CD [mack1998]
      vnlocal(6,1) = 'CT3'
      vnlocal(6,2) = 'CD'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.500d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT1 - CT1 [mack1998]
      vnlocal(1,1) = 'CT1'
      vnlocal(1,2) = 'CT1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CT2 - CA,CPB
c     --- CT3 - CA,CPB
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4900d0
      vclocal(1) = (230.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 4
      vnlocal(1,1) = 'CT2'
      vnlocal(1,2) = 'CA'
      vnlocal(2,1) = 'CT2'
      vnlocal(2,2) = 'CPB'
      vnlocal(3,1) = 'CT3'
      vnlocal(3,2) = 'CA'
      vnlocal(4,1) = 'CT3'
      vnlocal(4,2) = 'CPB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CT2,CT3 - CPH1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.500d0
      vclocal(1) = (229.63d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
      vnlocal(1,1) = 'CT2'
      vnlocal(1,2) = 'CPH1'
      vnlocal(2,1) = 'CT3'
      vnlocal(2,2) = 'CPH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5380d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CT2 - CT1 [mack1998]
      vnlocal(1,1) = 'CT2'
      vnlocal(1,2) = 'CT1'
c     --- CT3 - CT1 [mack1998]
      vnlocal(2,1) = 'CT3'
      vnlocal(2,2) = 'CT1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CT3 - CS
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5310d0
      vclocal(1) = (190.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CT3'
      vnlocal(1,2) = 'CS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5280d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT3 - CT2 [mack1998]
      vnlocal(1,1) = 'CT3'
      vnlocal(1,2) = 'CT2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5300d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CT3 - CT3 [mack1998]
      vnlocal(1,1) = 'CT3'
      vnlocal(1,2) = 'CT3'
c     --- CT2 - CT2 [mack1998]
      vnlocal(2,1) = 'CT2'
      vnlocal(2,2) = 'CT2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CY - CA
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (350.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CY'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CY - CPT
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4400d0
      vclocal(1) = (350.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CY'
      vnlocal(1,2) = 'CPT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CY - CT2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5100d0
      vclocal(1) = (230.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CY'
      vnlocal(1,2) = 'CT2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- FE - CM
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.900d0
      vclocal(1) = (258.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'FE'
      vnlocal(1,2) = 'CM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- FE - CPM
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 3.3814d0
      vclocal(1) = (0.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'FE'
      vnlocal(1,2) = 'CPM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- H - CD
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1100d0
      vclocal(1) = (330.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'H'
      vnlocal(1,2) = 'CD'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1000d0
      vclocal(1) = (330.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HA - C
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0830d0
      vclocal(1) = (340.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HA - CA
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- HA - CC
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1000d0
      vclocal(1) = (317.13d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'CC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1110d0
      vclocal(1) = (309.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 4
c     --- HA - CP2 [mack1998]
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'CP2'
c     --- HA - CP3 [mack1998]
      vnlocal(2,1) = 'HA'
      vnlocal(2,2) = 'CP3'
c     --- HA - CT1 [mack1998]
      vnlocal(3,1) = 'HA'
      vnlocal(3,2) = 'CT1'
c     --- HA - CT2 [mack1998]
      vnlocal(4,1) = 'HA'
      vnlocal(4,2) = 'CT2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- HA - CPM
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0900d0
      vclocal(1) = (367.60d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'CPM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- HA - CS
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1110d0
      vclocal(1) = (300.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'CS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1110d0
      vclocal(1) = (322.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- HA - CT3 [mack1998]
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'CT3'
c     --- CTL3 - HAL [mack1998]
      vnlocal(2,1) = 'CTL3'
      vnlocal(2,2) = 'HAL'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- HA - CY
c     --- HB - CP1,CT1,CT2,CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0800d0
      vclocal(1) = (330.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 5
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'CY'
      vnlocal(2,1) = 'HB'
      vnlocal(2,2) = 'CP1'
      vnlocal(3,1) = 'HB'
      vnlocal(3,2) = 'CT1'
      vnlocal(4,1) = 'HB'
      vnlocal(4,2) = 'CT2'
      vnlocal(5,1) = 'HB'
      vnlocal(5,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- HP - CA
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0800d0
      vclocal(1) = (340.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'HP'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0800d0
      vclocal(1) = (350.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HP - CY
      vnlocal(1,1) = 'HP'
      vnlocal(1,2) = 'CY'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- HR1 - CPH1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0830d0
      vclocal(1) = (375.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'HR1'
      vnlocal(1,2) = 'CPH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- HR1 - CPH2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0900d0
      vclocal(1) = (340.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'HR1'
      vnlocal(1,2) = 'CPH2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- HR2 - CPH2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0700d0
      vclocal(1) = (333.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'HR2'
      vnlocal(1,2) = 'CPH2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- HR3 - CPH1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0830d0
      vclocal(1) = (365.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'HR3'
      vnlocal(1,2) = 'CPH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5139d0
      vclocal(1) = (0.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HT - HT [mack1998] - only needed for shake, included here
c     ------------ for completeness
      vnlocal(1,1) = 'HT'
      vnlocal(1,2) = 'HT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- N - C
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3000d0
      vclocal(1) = (260.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'N'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- N - CP1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4340d0
      vclocal(1) = (320.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'N'
      vnlocal(1,2) = 'CP1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- N - CP3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4550d0
      vclocal(1) = (320.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'N'
      vnlocal(1,2) = 'CP3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NC2 - C
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3650d0
      vclocal(1) = (463.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NC2'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NC2 - CT2,CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4900d0
      vclocal(1) = (261.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
      vnlocal(1,1) = 'NC2'
      vnlocal(1,2) = 'CT2'
      vnlocal(2,1) = 'NC2'
      vnlocal(2,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NC2 - HC
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (455.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NC2'
      vnlocal(1,2) = 'HC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NH1 - C
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3450d0
      vclocal(1) = (370.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NH1'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NH1 - CT1,CT2,CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4300d0
      vclocal(1) = (320.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
      vnlocal(1,1) = 'NH1'
      vnlocal(1,2) = 'CT1'
      vnlocal(2,1) = 'NH1'
      vnlocal(2,2) = 'CT2'
      vnlocal(3,1) = 'NH1'
      vnlocal(3,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NH1 - H
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9970d0
      vclocal(1) = (440.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NH1'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9800d0
      vclocal(1) = (405.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- NH1 - HC
      vnlocal(1,1) = 'NH1'
      vnlocal(1,2) = 'HC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NH2 - CC
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3600d0
      vclocal(1) = (430.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NH2'
      vnlocal(1,2) = 'CC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NH2 - CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4550d0
      vclocal(1) = (240.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NH2'
      vnlocal(1,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NH2 - H
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (480.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NH2'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NH2 - HC
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (460.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NH2'
      vnlocal(1,2) = 'HC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NH3 - CT1,CT2,CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4800d0
      vclocal(1) = (200.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
      vnlocal(1,1) = 'NH3'
      vnlocal(1,2) = 'CT1'
      vnlocal(2,1) = 'NH3'
      vnlocal(2,2) = 'CT2'
      vnlocal(3,1) = 'NH3'
      vnlocal(3,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NH3 - HC
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0400d0
      vclocal(1) = (403.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NH3'
      vnlocal(1,2) = 'HC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NP - CP1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4850d0
      vclocal(1) = (320.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NP'
      vnlocal(1,2) = 'CP1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NP - CP3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5020d0
      vclocal(1) = (320.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NP'
      vnlocal(1,2) = 'CP3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NP - HC
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0060d0
      vclocal(1) = (460.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NP'
      vnlocal(1,2) = 'HC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NPH - CPA
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3757d0
      vclocal(1) = (377.20d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NPH'
      vnlocal(1,2) = 'CPA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NPH - FE
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.9580d0
      vclocal(1) = (270.20d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NPH'
      vnlocal(1,2) = 'FE'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NR1 - CPH1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3800d0
      vclocal(1) = (400.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NR1'
      vnlocal(1,2) = 'CPH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NR1 - CPH2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3600d0
      vclocal(1) = (400.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NR1'
      vnlocal(1,2) = 'CPH2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NR1 - H
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (466.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NR1'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NR2 - CPH1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3800d0
      vclocal(1) = (400.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NR2'
      vnlocal(1,2) = 'CPH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NR2 - CPH2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3200d0
      vclocal(1) = (400.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NR2'
      vnlocal(1,2) = 'CPH2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NR2 - FE
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.2000d0
      vclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NR2'
      vnlocal(1,2) = 'FE'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NR3 - CPH1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3700d0
      vclocal(1) = (380.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NR3'
      vnlocal(1,2) = 'CPH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NR3 - CPH2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3200d0
      vclocal(1) = (380.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NR3'
      vnlocal(1,2) = 'CPH2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NR3 - H
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (453.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NR3'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NY - CA
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3700d0
      vclocal(1) = (270.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NY'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NY - CPT
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3750d0
      vclocal(1) = (270.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NY'
      vnlocal(1,2) = 'CPT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- NY - H
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9760d0
      vclocal(1) = (465.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'NY'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- O - C
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.2300d0
      vclocal(1) = (620.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'O'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- O - CC
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.2300d0
      vclocal(1) = (650.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'O'
      vnlocal(1,2) = 'CC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.2200d0
      vclocal(1) = (750.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- OB - CC [mack1998]
      vnlocal(1,1) = 'OB'
      vnlocal(1,2) = 'CC'
c     --- OB - CD [mack1998]
      vnlocal(2,1) = 'OB'
      vnlocal(2,2) = 'CD'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.2600d0
      vclocal(1) = (525.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- OC - CA [mack1998]
      vnlocal(1,1) = 'OC'
      vnlocal(1,2) = 'CA'
c     --- OC - CC [mack1998]
      vnlocal(2,1) = 'OC'
      vnlocal(2,2) = 'CC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- OC - CT2,CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3300d0
      vclocal(1) = (450.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
      vnlocal(1,1) = 'OC'
      vnlocal(1,2) = 'CT2'
      vnlocal(2,1) = 'OC'
      vnlocal(2,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- OH1 - CA
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4110d0
      vclocal(1) = (334.30d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'OH1'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4000d0
      vclocal(1) = (230.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OH1 - CD [mack1998]
      vnlocal(1,1) = 'OH1'
      vnlocal(1,2) = 'CD'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4200d0
      vclocal(1) = (428.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- OH1 - CT1 [mack1998]
      vnlocal(1,1) = 'OH1'
      vnlocal(1,2) = 'CT1'
c     --- OH1 - CT2 [mack1998]
      vnlocal(2,1) = 'OH1'
      vnlocal(2,2) = 'CT2'
c     --- OH1 - CT3 [mack1998]
      vnlocal(3,1) = 'OH1'
      vnlocal(3,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9600d0
      vclocal(1) = (545.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OH1 - H [mack1998]
      vnlocal(1,1) = 'OH1'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1280d0
      vclocal(1) = (1115.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OM - CM
      vnlocal(1,1) = 'OM'
      vnlocal(1,2) = 'CM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- OM - FE
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.8000d0
      vclocal(1) = (250.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'OM'
      vnlocal(1,2) = 'FE'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- OM - OM
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.2300d0
      vclocal(1) = (600.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'OM'
      vnlocal(1,2) = 'OM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3340d0
      vclocal(1) = (150.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OS - CD [mack1998]
      vnlocal(1,1) = 'OS'
      vnlocal(1,2) = 'CD'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4300d0
      vclocal(1) = (340.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OS - CT3 [mack1998]
      vnlocal(1,1) = 'OS'
      vnlocal(1,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9572d0
      vclocal(1) = (450.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OT - HT [mack1998]
      vnlocal(1,1) = 'OT'
      vnlocal(1,2) = 'HT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.8180d0
      vclocal(1) = (198.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- S - CT2 [mack1998]
      vnlocal(1,1) = 'S'
      vnlocal(1,2) = 'CT2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- S - CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.8160d0
      vclocal(1) = (240.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'S'
      vnlocal(1,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- S - HS
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3250d0
      vclocal(1) = (275.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'S'
      vnlocal(1,2) = 'HS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- SM - CT2,CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.8160d0
      vclocal(1) = (214.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
      vnlocal(1,1) = 'SM'
      vnlocal(1,2) = 'CT2'
      vnlocal(2,1) = 'SM'
      vnlocal(2,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- SM - SM
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.0290d0
      vclocal(1) = (173.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'SM'
      vnlocal(1,2) = 'SM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- SM - CS
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.8360d0
      vclocal(1) = (205.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'SM'
      vnlocal(1,2) = 'CS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * charmm22 angles with urey-bradley                              *
c     * also includes many charmm27 retained from charmm22             *
c     ******************************************************************
c     --- MacKerell et al.
c     --- J. Phys. Chem. B (1998) 102, 3586-3616 supplemental info

      type = type + 1
c     --- CA - CA - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.0d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.4162d0
      bclocal(3) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CP1,CP3 - N - C
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 117.0d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CP1'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'C'
      bnlocal(2,1) = 'CP3'
      bnlocal(2,2) = 'N'
      bnlocal(2,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CP2 - CP1 - C,CC
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.30d0
      bclocal(1) = (52.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CP2'
      bnlocal(1,2) = 'CP1'
      bnlocal(1,3) = 'C'
      bnlocal(2,1) = 'CP2'
      bnlocal(2,2) = 'CP1'
      bnlocal(2,3) = 'CC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CP2 - CP1 - CD
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.30d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CP2'
      bnlocal(1,2) = 'CP1'
      bnlocal(1,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CP2 - CP2 - CP1
C     --- CP3 - CP2 - CP2
C     --- NP - CP1 - CP2
C     --- NP - CP3 - CP2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
      bnlocal(1,1) = 'CP2'
      bnlocal(1,2) = 'CP2'
      bnlocal(1,3) = 'CP1'
      bnlocal(2,1) = 'CP3'
      bnlocal(2,2) = 'CP2'
      bnlocal(2,3) = 'CP2'
      bnlocal(3,1) = 'NP'
      bnlocal(3,2) = 'CP1'
      bnlocal(3,3) = 'CP2'
      bnlocal(4,1) = 'NP'
      bnlocal(4,2) = 'CP3'
      bnlocal(4,3) = 'CP2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CP3 - N - CP1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 114.20d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CP3'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'CP1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CP3 - NP - CP1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.00d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CP3'
      bnlocal(1,2) = 'NP'
      bnlocal(1,3) = 'CP1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPA - CPB - C
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.74d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPA'
      bnlocal(1,2) = 'CPB'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPA - CPM - CPA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.12d0
      bclocal(1) = (94.20d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPA'
      bnlocal(1,2) = 'CPM'
      bnlocal(1,3) = 'CPA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPA - NPH - CPA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 103.90d0
      bclocal(1) = (139.30d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPA'
      bnlocal(1,2) = 'NPH'
      bnlocal(1,3) = 'CPA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 121.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CPB - C - C
      bnlocal(1,1) = 'CPB'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPB - CPB - C
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.75d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPB'
      bnlocal(1,2) = 'CPB'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPB - CPB - CPA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 106.51d0
      bclocal(1) = (30.80d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPB'
      bnlocal(1,2) = 'CPB'
      bnlocal(1,3) = 'CPA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPH2 - NR1 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.50d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPH2'
      bnlocal(1,2) = 'NR1'
      bnlocal(1,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPH2 - NR2 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 104.00d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPH2'
      bnlocal(1,2) = 'NR2'
      bnlocal(1,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPH2 - NR3 - CPH1
C     --- NR3 - CPH1 - CPH1
C     --- NR3 - CPH2 - NR3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.00d0
      bclocal(1) = (145.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'CPH2'
      bnlocal(1,2) = 'NR3'
      bnlocal(1,3) = 'CPH1'
      bnlocal(2,1) = 'NR3'
      bnlocal(2,2) = 'CPH1'
      bnlocal(2,3) = 'CPH1'
      bnlocal(3,1) = 'NR3'
      bnlocal(3,2) = 'CPH2'
      bnlocal(3,3) = 'NR3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPM - CPA - CPB
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 124.07d0
      bclocal(1) = (61.60d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPM'
      bnlocal(1,2) = 'CPA'
      bnlocal(1,3) = 'CPB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPT - CA - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 118.00d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPT'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPT - CPT - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.00d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPT'
      bnlocal(1,2) = 'CPT'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPT - CY - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.40d0
      bclocal(1) = (120.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2610d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPT'
      bnlocal(1,2) = 'CY'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CPT - NY -CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.00d0
      bclocal(1) = (110.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CPT'
      bnlocal(1,2) = 'NY'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.00d0
      bclocal(1) = (52.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
C     --- CT1 - CT1 - C [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT1'
      bnlocal(ii,2) = 'CT1'
      bnlocal(ii,3) = 'C'
C     --- CT1 - CT1 - CC [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT1'
      bnlocal(ii,2) = 'CT1'
      bnlocal(ii,3) = 'CC'
C     --- CT1 - CT2 - CC [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT1'
      bnlocal(ii,2) = 'CT2'
      bnlocal(ii,3) = 'CC'
C     --- CT1 - CT2 - CD [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT1'
      bnlocal(ii,2) = 'CT2'
      bnlocal(ii,3) = 'CD'
C     --- CT2 - CT1 - C [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT2'
      bnlocal(ii,2) = 'CT1'
      bnlocal(ii,3) = 'C'
C     --- CT2 - CT1 - CC [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT2'
      bnlocal(ii,2) = 'CT1'
      bnlocal(ii,3) = 'CC'
C     --- CT2 - CT1 - CD [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT2'
      bnlocal(ii,2) = 'CT1'
      bnlocal(ii,3) = 'CD'
c     --- CT2 - CT2 - C [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT2'
      bnlocal(ii,2) = 'CT2'
      bnlocal(ii,3) = 'C'
c     --- CT2 - CT2 - CC [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT2'
      bnlocal(ii,2) = 'CT2'
      bnlocal(ii,3) = 'CC'
c     --- CT2 - CT2 - CD [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT2'
      bnlocal(ii,2) = 'CT2'
      bnlocal(ii,3) = 'CD'
C     --- CT3 - CT1 - C [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT3'
      bnlocal(ii,2) = 'CT1'
      bnlocal(ii,3) = 'C'
C     --- CT3 - CT1 - CC [mack1998]
      ii = ii + 1
      bnlocal(ii,1) = 'CT3'
      bnlocal(ii,2) = 'CT1'
      bnlocal(ii,3) = 'CC'
c     --- MGM added similar CT3 - CT2 - CC on 9-29-2000
      ii = ii + 1
      bnlocal(ii,1) = 'CT3'
      bnlocal(ii,2) = 'CT2'
      bnlocal(ii,3) = 'CC'
c     --- MGM added similar CT3 - CT2 - CD on 08-27-2004
      ii = ii + 1
      bnlocal(ii,1) = 'CT3'
      bnlocal(ii,2) = 'CT2'
      bnlocal(ii,3) = 'CD'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.00d0
      bclocal(1) = (53.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CT1 - CT1 - CT1
      bnlocal(1,1) = 'CT1'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'CT1'
C     --- CT2 - CT1 - CT1
      bnlocal(2,1) = 'CT2'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CT1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT1 - CT2 - CA
C     --- CT3 - CT2 - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.50d0
      bclocal(1) = (51.80d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CT1'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CA'
      bnlocal(2,1) = 'CT3'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT1 - CT2 - CPH1
C     --- CT3 - CT2 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 113.00d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CT1'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CPH1'
      bnlocal(2,1) = 'CT3'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 113.50d0
      bclocal(1) = (58.350d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (11.16d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- CT1 - CT2 - CT1
      bnlocal(1,1) = 'CT1'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CT1'
C     --- CT2 - CT2 - CT1
      bnlocal(2,1) = 'CT2'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'CT1'
C     --- CT3 - CT2 - CT1
      bnlocal(3,1) = 'CT3'
      bnlocal(3,2) = 'CT2'
      bnlocal(3,3) = 'CT1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT1 - NH1 - C
C     --- CT2 - NH1 - C
C     --- CT3 - NH1 - C
C     --- H - NH2 - CC
C     --- HA - C - CPB
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 5
      bnlocal(1,1) = 'CT1'
      bnlocal(1,2) = 'NH1'
      bnlocal(1,3) = 'C'
      bnlocal(2,1) = 'CT2'
      bnlocal(2,2) = 'NH1'
      bnlocal(2,3) = 'C'
      bnlocal(3,1) = 'CT3'
      bnlocal(3,2) = 'NH1'
      bnlocal(3,3) = 'C'
      bnlocal(4,1) = 'H'
      bnlocal(4,2) = 'NH2'
      bnlocal(4,3) = 'CC'
      bnlocal(5,1) = 'HA'
      bnlocal(5,2) = 'C'
      bnlocal(5,3) = 'CPB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT2 - CA - CA
C     --- CT3 - CA - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.30d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CA'
      bnlocal(2,1) = 'CT3'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT2 - CPB - CPA
C     --- CT3 - CPB - CPA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.74d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'CPB'
      bnlocal(1,3) = 'CPA'
      bnlocal(2,1) = 'CT3'
      bnlocal(2,2) = 'CPB'
      bnlocal(2,3) = 'CPA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT2 - CPB - CPB
C     --- CT3 - CPB - CPB
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.75d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'CPB'
      bnlocal(1,3) = 'CPB'
      bnlocal(2,1) = 'CT3'
      bnlocal(2,2) = 'CPB'
      bnlocal(2,3) = 'CPB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT2 - CPH1 - CPH1
C     --- CT3 - CPH1 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 130.00d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'CPH1'
      bnlocal(2,1) = 'CT3'
      bnlocal(2,2) = 'CPH1'
      bnlocal(2,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT2 - CT2 - CPB
C     --- C - CH2E - CH2E [neria1996]
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 113.00d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CPB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 113.60d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (11.16d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CT2 - CT2 - CT2
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CT2'
C     --- CT2 - CT3 - CT1  (not clear how you can get this bend)
      bnlocal(2,1) = 'CT2'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'CT1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT2 - CY - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 129.40d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'CY'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT2 - CY - CPT
C     --- NR1 - CPH1 - CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 124.00d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'CY'
      bnlocal(1,3) = 'CPT'
      bnlocal(2,1) = 'NR1'
      bnlocal(2,2) = 'CPH1'
      bnlocal(2,3) = 'CT2'
      bnlocal(3,1) = 'NR1'
      bnlocal(3,2) = 'CPH1'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT2 - NC2 - C
C     --- CT3 - NC2 - C
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (62.30d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'NC2'
      bnlocal(1,3) = 'C'
      bnlocal(2,1) = 'CT3'
      bnlocal(2,2) = 'NC2'
      bnlocal(2,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.60d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2651d0
      bclocal(3) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CT2 - OS - CD [mack1998]
      bnlocal(1,1) = 'CT2'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'CD'
C     --- CT3 - OS - CD [mack1998]
      bnlocal(2,1) = 'CT3'
      bnlocal(2,2) = 'OS'
      bnlocal(2,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.50d0
      bclocal(1) = (53.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT3 - CT1 - CT1
      bnlocal(1,1) = 'CT3'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'CT1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT3 - CT1 - CT2
C     --- CT3 - CT1 - CT3
C     --- CT3 - CT2 - CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 114.00d0
      bclocal(1) = (53.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'CT3'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'CT2'
      bnlocal(2,1) = 'CT3'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CT3'
      bnlocal(3,1) = 'CT3'
      bnlocal(3,2) = 'CT2'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 115.00d0
      bclocal(1) = (58.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT3 - CT2 - CT2
      bnlocal(1,1) = 'CT3'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CT3 - S - CT2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 95.00d0
      bclocal(1) = (34.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CT3'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CY - CPT - CA
C     --- NY - CPT - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 130.60d0
      bclocal(1) = (160.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CY'
      bnlocal(1,2) = 'CPT'
      bnlocal(1,3) = 'CA'
      bnlocal(2,1) = 'NY'
      bnlocal(2,2) = 'CPT'
      bnlocal(2,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CY - CPT - CPT
C     --- NY - CPT - CPT
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.40d0
      bclocal(1) = (110.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CY'
      bnlocal(1,2) = 'CPT'
      bnlocal(1,3) = 'CPT'
      bnlocal(2,1) = 'NY'
      bnlocal(2,2) = 'CPT'
      bnlocal(2,3) = 'CPT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- CY - CT2 - CT1
C     --- CY - CT2 - CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 114.00d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CY'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CT1'
      bnlocal(2,1) = 'CY'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- FE - NPH - CPA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 128.05d0
      bclocal(1) = (96.15d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'FE'
      bnlocal(1,2) = 'NPH'
      bnlocal(1,3) = 'CPA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- FE - NR2 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 133.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'FE'
      bnlocal(1,2) = 'NR2'
      bnlocal(1,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- FE - NR2 - CPH2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 123.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'FE'
      bnlocal(1,2) = 'NR2'
      bnlocal(1,3) = 'CPH2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- H - NH1 - C
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 123.00d0
      bclocal(1) = (34.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NH1'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- H - NH1 - CT1,CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 117.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NH1'
      bnlocal(1,3) = 'CT1'
      bnlocal(2,1) = 'H'
      bnlocal(2,2) = 'NH1'
      bnlocal(2,3) = 'CT2'
      bnlocal(3,1) = 'H'
      bnlocal(3,2) = 'NH1'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- H - NH2 - H
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (23.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NH2'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- H - NR1 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.50d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1500d0
      bclocal(3) = (20.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NR1'
      bnlocal(1,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- H - NR1 - CPH2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 127.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1400d0
      bclocal(3) = (20.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NR1'
      bnlocal(1,3) = 'CPH2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- H - NR3 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1300d0
      bclocal(3) = (15.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NR3'
      bnlocal(1,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- H - NR3 - CPH2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.0900d0
      bclocal(3) = (15.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NR3'
      bnlocal(1,3) = 'CPH2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- H - NY - CA,CPT
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.00d0
      bclocal(1) = (28.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NY'
      bnlocal(1,3) = 'CA'
      bnlocal(2,1) = 'H'
      bnlocal(2,2) = 'NY'
      bnlocal(2,3) = 'CPT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- H - OH1 - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.00d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'OH1'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 115.00d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- H - OH1 - CD [mack1998]
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'OH1'
      bnlocal(1,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 106.00d0
      bclocal(1) = (57.50d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- H - OH1 - CT1
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'OH1'
      bnlocal(1,3) = 'CT1'
C     --- H - OH1 - CT2
      bnlocal(2,1) = 'H'
      bnlocal(2,2) = 'OH1'
      bnlocal(2,3) = 'CT2'
C     --- H - OH1 - CT3
      bnlocal(3,1) = 'H'
      bnlocal(3,2) = 'OH1'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - C - C
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.50d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - C - HA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 118.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - CA - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (29.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1525d0
      bclocal(3) = (25.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - CA - CPT
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.00d0
      bclocal(1) = (41.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CPT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - CA - CY
C     --- HP - CA - CY
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.00d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1730d0
      bclocal(3) = (25.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CY'
      bnlocal(2,1) = 'HP'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'CY'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (33.430d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- HA - CP2 - CP1 [mack1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CP2'
      bnlocal(1,3) = 'CP1'
C     --- HA - CT2 - CT1 [mack1998]
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'CT1'
C     --- HA - CT3 - CT1 [mack1998]
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CT3'
      bnlocal(3,3) = 'CT1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (26.50d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
C     --- HA - CP2 - CP2
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CP2'
      bnlocal(1,3) = 'CP2'
C     --- HA - CP2 - CP3
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CP2'
      bnlocal(2,3) = 'CP3'
C     --- HA - CP3 - CP2
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CP3'
      bnlocal(3,3) = 'CP2'
C     --- HA - CT2 - CT2
      bnlocal(4,1) = 'HA'
      bnlocal(4,2) = 'CT2'
      bnlocal(4,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.00d0
      bclocal(1) = (35.50d0)*(CON_KCALTOK)
      bclocal(2) = 1.8020d0
      bclocal(3) = (5.40d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
C     --- HA - CP2 - HA [mack1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CP2'
      bnlocal(1,3) = 'HA'
C     --- HA - CP3 - HA [mack1998]
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CP3'
      bnlocal(2,3) = 'HA'
C     --- HA - CT1 - HA [mack1998] (not sure how you can get this one)
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CT1'
      bnlocal(3,3) = 'HA'
C     --- HA - CT2 - HA [mack1998]
      bnlocal(4,1) = 'HA'
      bnlocal(4,2) = 'CT2'
      bnlocal(4,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - CPM - CPA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 117.44d0
      bclocal(1) = (12.70d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CPM'
      bnlocal(1,3) = 'CPA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - CPM - FE
C     --- OM - OM - FE
C     --- FE - OM - OM [neria1996]
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 180.00d0
      bclocal(1) = (0.0d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CPM'
      bnlocal(1,3) = 'FE'
      bnlocal(2,1) = 'OM'
      bnlocal(2,2) = 'OM'
      bnlocal(2,3) = 'FE'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (34.60d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
C     --- HA - CS - CT3
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CS'
      bnlocal(1,3) = 'CT3'
C     --- HA - CT2 - CT3
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'CT3'
C     --- HA - CT3 - CS
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CT3'
      bnlocal(3,3) = 'CS'
C     --- HA - CT3 - CT2
      bnlocal(4,1) = 'HA'
      bnlocal(4,2) = 'CT3'
      bnlocal(4,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - CS - HA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.40d0
      bclocal(1) = (35.50d0)*(CON_KCALTOK)
      bclocal(2) = 1.7750d0
      bclocal(3) = (14.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CS'
      bnlocal(1,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (33.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1630d0
      bclocal(3) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 8
C     --- HA - CT1 - C [mack1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'C'
C     --- HA - CT1 - CD [mack1998]
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CD'
C     --- HA - CT2 - C [mack1998]
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CT2'
      bnlocal(3,3) = 'C'
C     --- HA - CT2 - CC [mack1998]
      bnlocal(4,1) = 'HA'
      bnlocal(4,2) = 'CT2'
      bnlocal(4,3) = 'CC'
C     --- HA - CT2 - CD [mack1998]
      bnlocal(5,1) = 'HA'
      bnlocal(5,2) = 'CT2'
      bnlocal(5,3) = 'CD'
C     --- HA - CT3 - C [mack1998]
      bnlocal(6,1) = 'HA'
      bnlocal(6,2) = 'CT3'
      bnlocal(6,3) = 'C'
C     --- HA - CT3 - CC [mack1998]
      bnlocal(7,1) = 'HA'
      bnlocal(7,2) = 'CT3'
      bnlocal(7,3) = 'CC'
C     --- HA - CT3 - CD [mack1998]
      bnlocal(8,1) = 'HA'
      bnlocal(8,2) = 'CT3'
      bnlocal(8,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (34.50d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- HA - CT1 - CT1
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'CT1'
C     --- HA - CT1 - CT2
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CT2'
C     --- HA - CT1 - CT3
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CT1'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.50d0
      bclocal(1) = (49.30d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- HA - CT2 - CA
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CA'
C     --- HA - CT3 - CA
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 9
C     --- HA - CT2 - CPB
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CPB'
C     --- HA - CT3 - CPB
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'CPB'
C     --- HB - CT1 - C
      bnlocal(3,1) = 'HB'
      bnlocal(3,2) = 'CT1'
      bnlocal(3,3) = 'C'
C     --- HB - CT1 - CC
      bnlocal(4,1) = 'HB'
      bnlocal(4,2) = 'CT1'
      bnlocal(4,3) = 'CC'
C     --- HB - CT1 - CD
      bnlocal(5,1) = 'HB'
      bnlocal(5,2) = 'CT1'
      bnlocal(5,3) = 'CD'
C     --- HB - CT2 - C
      bnlocal(6,1) = 'HB'
      bnlocal(6,2) = 'CT2'
      bnlocal(6,3) = 'C'
C     --- HB - CT2 - CC
      bnlocal(7,1) = 'HB'
      bnlocal(7,2) = 'CT2'
      bnlocal(7,3) = 'CC'
C     --- HB - CT2 - CD
      bnlocal(8,1) = 'HB'
      bnlocal(8,2) = 'CT2'
      bnlocal(8,3) = 'CD'
C     --- HB - CT3 - C
      bnlocal(9,1) = 'HB'
      bnlocal(9,2) = 'CT3'
      bnlocal(9,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - CT2 - CPH1,CY
C     --- HA - CT3 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (33.43d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CPH1'
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'CY'
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CT3'
      bnlocal(3,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (37.50d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CT3 - CT3
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CT3'
      bnlocal(1,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.40d0
      bclocal(1) = (35.50d0)*(CON_KCALTOK)
      bclocal(2) = 1.8020d0
      bclocal(3) = (5.40d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CT3 - HA [mack1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CT3'
      bnlocal(1,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - CY - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.40d0
      bclocal(1) = (20.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1860d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CY'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HA - CY - CPT
C     --- HP - CY - CPT
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.40d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2550d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CY'
      bnlocal(1,3) = 'CPT'
      bnlocal(2,1) = 'HP'
      bnlocal(2,2) = 'CY'
      bnlocal(2,3) = 'CPT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HB - CP1 - C,CC,CD
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'HB'
      bnlocal(1,2) = 'CP1'
      bnlocal(1,3) = 'C'
      bnlocal(2,1) = 'HB'
      bnlocal(2,2) = 'CP1'
      bnlocal(2,3) = 'CC'
      bnlocal(3,1) = 'HB'
      bnlocal(3,2) = 'CP1'
      bnlocal(3,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HB - CP1 - CP2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 118.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HB'
      bnlocal(1,2) = 'CP1'
      bnlocal(1,3) = 'CP2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HB - CT1 - CT1,CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'HB'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'CT1'
      bnlocal(2,1) = 'HB'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CT2'
      bnlocal(3,1) = 'HB'
      bnlocal(3,2) = 'CT1'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HB - CT2 - HB
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 115.00d0
      bclocal(1) = (36.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HB'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'HB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HC - NC2 - C
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (49.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NC2'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HC - NC2 - CT2,CT3
C     --- OC - CA - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (40.40d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NC2'
      bnlocal(1,3) = 'CT2'
      bnlocal(2,1) = 'HC'
      bnlocal(2,2) = 'NC2'
      bnlocal(2,3) = 'CT3'
      bnlocal(3,1) = 'OC'
      bnlocal(3,2) = 'CA'
      bnlocal(3,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HC - NC2 - HC
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NC2'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HC - NH2 - CT3
C     --- CH3E - CH1E - CH3E [neria1996]
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NH2'
      bnlocal(1,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HC - NH2 - HC
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 106.50d0
      bclocal(1) = (39.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NH2'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HC - NH3 - CT1,CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.0740d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NH3'
      bnlocal(1,3) = 'CT1'
      bnlocal(2,1) = 'HC'
      bnlocal(2,2) = 'NH3'
      bnlocal(2,3) = 'CT2'
      bnlocal(3,1) = 'HC'
      bnlocal(3,2) = 'NH3'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HC - NH3 - HC
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NH3'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (33.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.0560d0
      bclocal(3) = (4.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- HC - NP - CP1
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NP'
      bnlocal(1,3) = 'CP1'
C     --- HC - NP - CP3
      bnlocal(2,1) = 'HC'
      bnlocal(2,2) = 'NP'
      bnlocal(2,3) = 'CP3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HC - NP - HC
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.50d0
      bclocal(1) = (51.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NP'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HP - CA - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1525d0
      bclocal(3) = (22.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HP'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HP - CA - CPT
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1460d0
      bclocal(3) = (22.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HP'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CPT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HP - CY - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.40d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1860d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HP'
      bnlocal(1,2) = 'CY'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HR1 - CPH1 - CPH1 
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 130.00d0
      bclocal(1) = (22.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2150d0
      bclocal(3) = (15.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HR1'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HR3 - CPH1 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 130.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.200d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HR3'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HS - S - CT2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 95.00d0
      bclocal(1) = (38.80d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HS'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- HS - S - CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 95.00d0
      bclocal(1) = (43.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HS'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 104.52d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HT - OT - HT
      bnlocal(1,1) = 'HT'
      bnlocal(1,2) = 'OT'
      bnlocal(1,3) = 'HT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- N - C - CP1,CT1,CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.50d0
      bclocal(1) = (20.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CP1'
      bnlocal(2,1) = 'N'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'CT1'
      bnlocal(3,1) = 'N'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'CT2'
      bnlocal(4,1) = 'N'
      bnlocal(4,2) = 'C'
      bnlocal(4,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- N - CP1 - C,CC,CD
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.20d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'CP1'
      bnlocal(1,3) = 'C'
      bnlocal(2,1) = 'N'
      bnlocal(2,2) = 'CP1'
      bnlocal(2,3) = 'CC'
      bnlocal(3,1) = 'N'
      bnlocal(3,2) = 'CP1'
      bnlocal(3,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- N - CP1 - CP2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.80d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'CP1'
      bnlocal(1,3) = 'CP2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- N - CP1 - HB
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.00d0
      bclocal(1) = (48.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'CP1'
      bnlocal(1,3) = 'HB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- N - CP3 - CP2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'CP3'
      bnlocal(1,3) = 'CP2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- N - CP3 - HA
C     --- NH1 - CT1 - HB
C     --- NH1 - CT2 - HB
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.00d0
      bclocal(1) = (48.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'CP3'
      bnlocal(1,3) = 'HA'
      bnlocal(2,1) = 'NH1'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'HB'
      bnlocal(3,1) = 'NH1'
      bnlocal(3,2) = 'CT2'
      bnlocal(3,3) = 'HB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NC2 - C - NC2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (52.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.3642d0
      bclocal(3) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NC2'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NC2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NC2 - CT2 - CT2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.50d0
      bclocal(1) = (67.70d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NC2'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NC2 - CT2 - HA
C     --- NC2 - CT3 - HA
C     --- NH3 - CT1 - HB
C     --- NH3 - CT2 - HB
C     --- NP - CP1 - HB
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.50d0
      bclocal(1) = (51.50d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 5
      bnlocal(1,1) = 'NC2'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'HA'
      bnlocal(2,1) = 'NC2'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'HA'
      bnlocal(3,1) = 'NH3'
      bnlocal(3,2) = 'CT1'
      bnlocal(3,3) = 'HB'
      bnlocal(4,1) = 'NH3'
      bnlocal(4,2) = 'CT2'
      bnlocal(4,3) = 'HB'
      bnlocal(5,1) = 'NP'
      bnlocal(5,2) = 'CP1'
      bnlocal(5,3) = 'HB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NH1 - C - CP1,CT1,CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 116.50d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
      bnlocal(1,1) = 'NH1'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CP1'
      bnlocal(2,1) = 'NH1'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'CT1'
      bnlocal(3,1) = 'NH1'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'CT2'
      bnlocal(4,1) = 'NH1'
      bnlocal(4,2) = 'C'
      bnlocal(4,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NH1 - CT1 - C,CC,CD
C     --- NH1 - CT2 - C,CC,CD
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 6
      bnlocal(1,1) = 'NH1'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'C'
      bnlocal(2,1) = 'NH1'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CC'
      bnlocal(3,1) = 'NH1'
      bnlocal(3,2) = 'CT1'
      bnlocal(3,3) = 'CD'
      bnlocal(4,1) = 'NH1'
      bnlocal(4,2) = 'CT2'
      bnlocal(4,3) = 'C'
      bnlocal(5,1) = 'NH1'
      bnlocal(5,2) = 'CT2'
      bnlocal(5,3) = 'CC'
      bnlocal(6,1) = 'NH1'
      bnlocal(6,2) = 'CT2'
      bnlocal(6,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 113.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
C     --- NH1 - CT1 - CT1
      bnlocal(1,1) = 'NH1'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'CT1'
C     --- NH1 - CT1 - CT2
      bnlocal(2,1) = 'NH1'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CT2'
C     --- NH1 - CT1 - CT3
      bnlocal(3,1) = 'NH1'
      bnlocal(3,2) = 'CT1'
      bnlocal(3,3) = 'CT3'
C     --- NH1 - CT2 - CT2
      bnlocal(4,1) = 'NH1'
      bnlocal(4,2) = 'CT2'
      bnlocal(4,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (51.50d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- NH1 - CT2 - HA [mack1998]
      bnlocal(1,1) = 'NH1'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'HA'
C     --- NH1 - CT3 - HA [mack1998]
      bnlocal(2,1) = 'NH1'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NH2 - CC - CP1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.50d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NH2'
      bnlocal(1,2) = 'CC'
      bnlocal(1,3) = 'CP1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NH2 - CC - CT1,CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 116.50d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.4500d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'NH2'
      bnlocal(1,2) = 'CC'
      bnlocal(1,3) = 'CT1'
      bnlocal(2,1) = 'NH2'
      bnlocal(2,2) = 'CC'
      bnlocal(2,3) = 'CT2'
      bnlocal(3,1) = 'NH2'
      bnlocal(3,2) = 'CC'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.00d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      bclocal(2) = 1.9800d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- NH2 - CC - HA
      bnlocal(1,1) = 'NH2'
      bnlocal(1,2) = 'CC'
      bnlocal(1,3) = 'HA'
C     --- NH1 - C - HA (CHO parameters for formamide cap on gramicidin)
c     --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c         `High-resolution polypeptide structure in a lamellar 
c          phase lipid environment from solid state NMR derived 
c          orientational constraints'
c         Structure, 5, 1655-1669 (1997).
c     --- T.B. Woolf personal communication 08-16-2001 (M.G. Martin)
      bnlocal(2,1) = 'NH1'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1400d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NH2 - CT3 - HA
      bnlocal(1,1) = 'NH2'
      bnlocal(1,2) = 'CT3'
      bnlocal(1,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NH3 - CT1 - C,CC
C     --- NH3 - CT2 - C,CC,CD
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.00d0
      bclocal(1) = (43.70d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 5
      bnlocal(1,1) = 'NH3'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'C'
      bnlocal(2,1) = 'NH3'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CC'
      bnlocal(3,1) = 'NH3'
      bnlocal(3,2) = 'CT2'
      bnlocal(3,3) = 'C'
      bnlocal(4,1) = 'NH3'
      bnlocal(4,2) = 'CT2'
      bnlocal(4,3) = 'CC'
      bnlocal(5,1) = 'NH3'
      bnlocal(5,2) = 'CT2'
      bnlocal(5,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.00d0
      bclocal(1) = (67.70d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
C     --- NH3 - CT1 - CT1
      bnlocal(1,1) = 'NH3'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'CT1'
C     --- NH3 - CT1 - CT2
      bnlocal(2,1) = 'NH3'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CT2'
C     --- NH3 - CT1 - CT3
      bnlocal(3,1) = 'NH3'
      bnlocal(3,2) = 'CT1'
      bnlocal(3,3) = 'CT3'
C     --- NH3 - CT2 - CT2
      bnlocal(4,1) = 'NH3'
      bnlocal(4,2) = 'CT2'
      bnlocal(4,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NH3 - CT2 - HA
C     --- NH3 - CT3 - HA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.50d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1010d0
      bclocal(3) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'NH3'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'HA'
      bnlocal(2,1) = 'NH3'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NP - CP1 - C,CC,CD
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 106.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'NP'
      bnlocal(1,2) = 'CP1'
      bnlocal(1,3) = 'C'
      bnlocal(2,1) = 'NP'
      bnlocal(2,2) = 'CP1'
      bnlocal(2,3) = 'CC'
      bnlocal(3,1) = 'NP'
      bnlocal(3,2) = 'CP1'
      bnlocal(3,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NP - CP3 - HA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.15d0
      bclocal(1) = (51.50d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NP'
      bnlocal(1,2) = 'CP3'
      bnlocal(1,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NPH - CPA - CPB
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.54d0
      bclocal(1) = (122.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NPH'
      bnlocal(1,2) = 'CPA'
      bnlocal(1,3) = 'CPB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NPH - CPA - CPM
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 124.39d0
      bclocal(1) = (88.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NPH'
      bnlocal(1,2) = 'CPA'
      bnlocal(1,3) = 'CPM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NPH - FE - CM
C     --- NR2 - FE - NPH
C     --- NP - FE - NP [neria1996]
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 90.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'NPH'
      bnlocal(1,2) = 'FE'
      bnlocal(1,3) = 'CM'
      bnlocal(2,1) = 'NR2'
      bnlocal(2,2) = 'FE'
      bnlocal(2,3) = 'NPH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NPH - FE - CPM
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 45.00d0
      bclocal(1) = (0.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NPH'
      bnlocal(1,2) = 'FE'
      bnlocal(1,3) = 'CPM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NPH - FE - NPH
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 90.00d0
      bclocal(1) = (14.39d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NPH'
      bnlocal(1,2) = 'FE'
      bnlocal(1,3) = 'NPH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR1 - CPH1 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 106.00d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR1'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR1 - CPH1 - HR3,HR1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 124.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1400d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'NR1'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'HR3'
      bnlocal(2,1) = 'NR1'
      bnlocal(2,2) = 'CPH1'
      bnlocal(2,3) = 'HR1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR2 - CPH1 - CPH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.00d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR2'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'CPH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR2 - CPH1 - CT2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR2'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR2 - CPH1 - HR3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1400d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR2'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'HR3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR2 - CPH2 - HR1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1200d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR2'
      bnlocal(1,2) = 'CPH2'
      bnlocal(1,3) = 'HR1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR2 - CPH2 - NR1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.50d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR2'
      bnlocal(1,2) = 'CPH2'
      bnlocal(1,3) = 'NR1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR2 - FE - CM
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 180.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR2'
      bnlocal(1,2) = 'FE'
      bnlocal(1,3) = 'CM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR3 - CPH1 - CT2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.00d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR3'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'CT2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR3 - CPH1 - HR1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.00d0
      bclocal(1) = (22.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1800d0
      bclocal(3) = (15.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR3'
      bnlocal(1,2) = 'CPH1'
      bnlocal(1,3) = 'HR1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR3 - CPH2 - HR2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.00d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1400d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR3'
      bnlocal(1,2) = 'CPH2'
      bnlocal(1,3) = 'HR2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NY - CA - CY
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.00d0
      bclocal(1) = (120.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2400d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NY'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CY'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NY - CA - HA,HP
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.00d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1770d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'NY'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'HA'
      bnlocal(2,1) = 'NY'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'HP'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- O - C - CP1
C     --- O - CC - CP1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 118.00d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CP1'
      bnlocal(2,1) = 'O'
      bnlocal(2,2) = 'CC'
      bnlocal(2,3) = 'CP1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- O - C - CT1,CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 121.00d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CT1'
      bnlocal(2,1) = 'O'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'CT2'
      bnlocal(3,1) = 'O'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 121.70d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- O - C - H
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'H'
C     --- O - C - HA (ethanolamine cap) References for this group are
c     --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c         `High-resolution polypeptide structure in a lamellar 
c          phase lipid environment from solid state NMR derived 
c          orientational constraints'
c         Structure, 5, 1655-1669 (1997).
c     --- T.B. Woolf personal communication 08-16-2001 (M.G. Martin)
      bnlocal(2,1) = 'O'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- O - C - N,NH1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.50d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'N'
      bnlocal(2,1) = 'O'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'NH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- O - CC - CT1,CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 121.00d0
      bclocal(1) = (15.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.4400d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'CC'
      bnlocal(1,3) = 'CT1'
      bnlocal(2,1) = 'O'
      bnlocal(2,2) = 'CC'
      bnlocal(2,3) = 'CT2'
      bnlocal(3,1) = 'O'
      bnlocal(3,2) = 'CC'
      bnlocal(3,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- O - CC - HA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.00d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'CC'
      bnlocal(1,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- O - CC - NH2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.50d0
      bclocal(1) = (75.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.3700d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'CC'
      bnlocal(1,3) = 'NH2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.00d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.4420d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
C     --- OB - CD - CP1 [mack1998]
      bnlocal(1,1) = 'OB'
      bnlocal(1,2) = 'CD'
      bnlocal(1,3) = 'CP1'
C     --- OB - CD - CT1 [mack1998]
      bnlocal(2,1) = 'OB'
      bnlocal(2,2) = 'CD'
      bnlocal(2,3) = 'CT1'
C     --- OB - CD - CT2 [mack1998]
      bnlocal(3,1) = 'OB'
      bnlocal(3,2) = 'CD'
      bnlocal(3,3) = 'CT2'
C     --- OB - CD - CT3 [mack1998]
      bnlocal(4,1) = 'OB'
      bnlocal(4,2) = 'CD'
      bnlocal(4,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 118.00d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.3880d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
C     --- OC - CC - CP1 [mack1998]
      bnlocal(1,1) = 'OC'
      bnlocal(1,2) = 'CC'
      bnlocal(1,3) = 'CP1'
C     --- OC - CC - CT1 [mack1998]
      bnlocal(2,1) = 'OC'
      bnlocal(2,2) = 'CC'
      bnlocal(2,3) = 'CT1'
C     --- OC - CC - CT2 [mack1998]
      bnlocal(3,1) = 'OC'
      bnlocal(3,2) = 'CC'
      bnlocal(3,3) = 'CT2'
C     --- OC - CC - CT3 [mack1998]
      bnlocal(4,1) = 'OC'
      bnlocal(4,2) = 'CC'
      bnlocal(4,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 124.00d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2250d0
      bclocal(3) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- OC - CC - OC
      bnlocal(1,1) = 'OC'
      bnlocal(1,2) = 'CC'
      bnlocal(1,3) = 'OC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- OC - CT2 - CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.00d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'OC'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- OC - CT2 - HA
C     --- OC - CT3 - HA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 118.30d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'OC'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'HA'
      bnlocal(2,1) = 'OC'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- OH1 - CA - CA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (45.20d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'OH1'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.50d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- OH1 - CD - CT2 [mack1998]
      bnlocal(1,1) = 'OH1'
      bnlocal(1,2) = 'CD'
      bnlocal(1,3) = 'CT2'
C     --- OH1 - CD - CT3 [mack1998]
      bnlocal(2,1) = 'OH1'
      bnlocal(2,2) = 'CD'
      bnlocal(2,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 123.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2620d0
      bclocal(3) = (210.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- OH1 - CD - OB
      bnlocal(1,1) = 'OH1'
      bnlocal(1,2) = 'CD'
      bnlocal(1,3) = 'OB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (75.70d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 5
C     --- OH1 - CT1 - CT1 [mack1998]
      bnlocal(1,1) = 'OH1'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'CT1'
C     --- OH1 - CT1 - CT3 [mack1998]
      bnlocal(2,1) = 'OH1'
      bnlocal(2,2) = 'CT1'
      bnlocal(2,3) = 'CT3'
C     --- OH1 - CT2 - CT1 [mack1998]
      bnlocal(3,1) = 'OH1'
      bnlocal(3,2) = 'CT2'
      bnlocal(3,3) = 'CT1'
C     --- OH1 - CT2 - CT2 [mack1998]
      bnlocal(4,1) = 'OH1'
      bnlocal(4,2) = 'CT2'
      bnlocal(4,3) = 'CT2'
C     --- OH1 - CT2 - CT3 [mack1998]
      bnlocal(5,1) = 'OH1'
      bnlocal(5,2) = 'CT2'
      bnlocal(5,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.89d0
      bclocal(1) = (45.90d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- OH1 - CT1 - HA
      bnlocal(1,1) = 'OH1'
      bnlocal(1,2) = 'CT1'
      bnlocal(1,3) = 'HA'
C     --- OH1 - CT2 - HA
      bnlocal(2,1) = 'OH1'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'HA'
C     --- OH1 - CT3 - HA
      bnlocal(3,1) = 'OH1'
      bnlocal(3,2) = 'CT3'
      bnlocal(3,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- OM - CM - FE
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 180.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'OM'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'FE'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- OM - FE - NPH
C     --- FE - CM - OM [neria1996]
C     --- FE - FE - CM [neria1996]
C     --- NP - FE - 0M [neria1996]
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 90.00d0
      bclocal(1) = (5.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'OM'
      bnlocal(1,2) = 'FE'
      bnlocal(1,3) = 'NPH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.00d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.3260d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
C     --- OS - CD - CP1 [mack1998]
      bnlocal(1,1) = 'OS'
      bnlocal(1,2) = 'CD'
      bnlocal(1,3) = 'CP1'
C     --- OS - CD - CT1 [mack1998]
      bnlocal(2,1) = 'OS'
      bnlocal(2,2) = 'CD'
      bnlocal(2,3) = 'CT1'
C     --- OS - CD - CT2 [mack1998]
      bnlocal(3,1) = 'OS'
      bnlocal(3,2) = 'CD'
      bnlocal(3,3) = 'CT2'
C     --- OS - CD - CT3 [mack1998]
      bnlocal(4,1) = 'OS'
      bnlocal(4,2) = 'CD'
      bnlocal(4,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.90d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2576d0
      bclocal(3) = (160.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- OS - CD - OB [mack1998]
      bnlocal(1,1) = 'OS'
      bnlocal(1,2) = 'CD'
      bnlocal(1,3) = 'OB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- OS - CT2 - HA
      bnlocal(1,1) = 'OS'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'HA'
C     --- OS - CT3 - HA
      bnlocal(2,1) = 'OS'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- S - CT2 - CT1
C     --- SM - CT2 - CT1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.50d0
      bclocal(1) = (58.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'S'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CT1'
      bnlocal(2,1) = 'SM'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'CT1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- S - CT2 - CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 114.50d0
      bclocal(1) = (58.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'S'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'CT2'
      bnlocal(2,1) = 'S'
      bnlocal(2,2) = 'CT2'
      bnlocal(2,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- S - CT2 - HA
C     --- S - CT3 - HA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.30d0
      bclocal(1) = (46.10d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'S'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'HA'
      bnlocal(2,1) = 'S'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- SM - CT2 - HA
C     --- SM - CT3 - HA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.00d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'SM'
      bnlocal(1,2) = 'CT2'
      bnlocal(1,3) = 'HA'
      bnlocal(2,1) = 'SM'
      bnlocal(2,2) = 'CT3'
      bnlocal(2,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- SM - SM - CT2,CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 103.30d0
      bclocal(1) = (72.50d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'SM'
      bnlocal(1,2) = 'SM'
      bnlocal(1,3) = 'CT2'
      bnlocal(2,1) = 'SM'
      bnlocal(2,2) = 'SM'
      bnlocal(2,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- SS - CS - CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 118.00d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'SS'
      bnlocal(1,2) = 'CS'
      bnlocal(1,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- SS - CS - HA
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.30d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'SS'
      bnlocal(1,2) = 'CS'
      bnlocal(1,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
C     --- NR1 - CPH2 - HR1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.50d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1400d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'NR1'
      bnlocal(1,2) = 'CPH2'
      bnlocal(1,3) = 'HR1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Charmm22 standard dihedral                                     *
c     * plus those charmm27 parameters that did not change from 22     *
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.200d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 6
c     --- C - CT1 - NH1 - C
      tnlocal(1,1) = 'C'
      tnlocal(1,2) = 'CT1'
      tnlocal(1,3) = 'NH1'
      tnlocal(1,4) = 'C'
c     --- C - CT2 - NH1 - C
      tnlocal(2,1) = 'C'
      tnlocal(2,2) = 'CT2'
      tnlocal(2,3) = 'NH1'
      tnlocal(2,4) = 'C'
c     --- CC - CT1 - NH1 - C
      tnlocal(3,1) = 'CC'
      tnlocal(3,2) = 'CT1'
      tnlocal(3,3) = 'NH1'
      tnlocal(3,4) = 'C'
c     --- CC - CT2 - NH1 - C
      tnlocal(4,1) = 'CC'
      tnlocal(4,2) = 'CT2'
      tnlocal(4,3) = 'NH1'
      tnlocal(4,4) = 'C'
c     --- CD - CT1 - NH1 - C
      tnlocal(5,1) = 'CD'
      tnlocal(5,2) = 'CT1'
      tnlocal(5,3) = 'NH1'
      tnlocal(5,4) = 'C'
c     --- CD - CT2 - NH1 - C
      tnlocal(6,1) = 'CD'
      tnlocal(6,2) = 'CT2'
      tnlocal(6,3) = 'NH1'
      tnlocal(6,4) = 'C'
      call twh_set_tornames(type,ntslocal,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.800d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- C - N - CP1 - C
      tnlocal(1,1) = 'C'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CP1'
      tnlocal(1,4) = 'C'
c     --- CC - CP1 - N - C
      tnlocal(2,1) = 'CC'
      tnlocal(2,2) = 'CP1'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'C'
c     --- CP2 - CP1 - N - C
      tnlocal(3,1) = 'CP2'
      tnlocal(3,2) = 'CP1'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'C'
c     --- HB - CP1 - N - C
      tnlocal(4,1) = 'HB'
      tnlocal(4,2) = 'CP1'
      tnlocal(4,3) = 'N'
      tnlocal(4,4) = 'C'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.100d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 8
c     --- CA - CA - CA - CA
      tnlocal(1,1) = 'CA'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'CA'
      tnlocal(1,4) = 'CA'
c     --- CA - CPT - CPT - CA
      tnlocal(2,1) = 'CA'
      tnlocal(2,2) = 'CPT'
      tnlocal(2,3) = 'CPT'
      tnlocal(2,4) = 'CA'
c     --- CPT - CA - CA - CA
      tnlocal(3,1) = 'CPT'
      tnlocal(3,2) = 'CA'
      tnlocal(3,3) = 'CA'
      tnlocal(3,4) = 'CA'
c     --- CPT - CPT - CA - CA
      tnlocal(4,1) = 'CPT'
      tnlocal(4,2) = 'CPT'
      tnlocal(4,3) = 'CA'
      tnlocal(4,4) = 'CA'
c     --- CT2 - CA - CA - CA
      tnlocal(5,1) = 'CT2'
      tnlocal(5,2) = 'CA'
      tnlocal(5,3) = 'CA'
      tnlocal(5,4) = 'CA'
c     --- CT3 - CA - CA - CA
      tnlocal(6,1) = 'CT3'
      tnlocal(6,2) = 'CA'
      tnlocal(6,3) = 'CA'
      tnlocal(6,4) = 'CA'
c     --- OC - CA - CA - CA
      tnlocal(7,1) = 'OC'
      tnlocal(7,2) = 'CA'
      tnlocal(7,3) = 'CA'
      tnlocal(7,4) = 'CA'
c     --- OH1 - CA - CA - CA
      tnlocal(8,1) = 'OH1'
      tnlocal(8,2) = 'CA'
      tnlocal(8,3) = 'CA'
      tnlocal(8,4) = 'CA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.040d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- CA - CT2 - CT1 - C
      tnlocal(1,1) = 'CA'
      tnlocal(1,2) = 'CT2'
      tnlocal(1,3) = 'CT1'
      tnlocal(1,4) = 'C'
c     --- CC - CT1 - CT2 - CA
      tnlocal(2,1) = 'CC'
      tnlocal(2,2) = 'CT1'
      tnlocal(2,3) = 'CT2'
      tnlocal(2,4) = 'CA'
c     --- HA - CT1 - CT2 - CA
      tnlocal(3,1) = 'HA'
      tnlocal(3,2) = 'CT1'
      tnlocal(3,3) = 'CT2'
      tnlocal(3,4) = 'CA'
c     --- HA - CT3 - CT2 - CA
      tnlocal(4,1) = 'HA'
      tnlocal(4,2) = 'CT3'
      tnlocal(4,3) = 'CT2'
      tnlocal(4,4) = 'CA'
c     --- N - CT1 - CT2 - CA
      tnlocal(5,1) = 'N'
      tnlocal(5,2) = 'CT1'
      tnlocal(5,3) = 'CT2'
      tnlocal(5,4) = 'CA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 31
c     --- CA - CY - CPT - CA
      tnlocal(1,1) = 'CA'
      tnlocal(1,2) = 'CY'
      tnlocal(1,3) = 'CPT'
      tnlocal(1,4) = 'CA'
c     --- CA - NY - CPT - CA
      tnlocal(2,1) = 'CA'
      tnlocal(2,2) = 'NY'
      tnlocal(2,3) = 'CPT'
      tnlocal(2,4) = 'CA'
c     --- CT2 - CPH1 - NR1 - CPH2
      tnlocal(3,1) = 'CT2'
      tnlocal(3,2) = 'CPH1'
      tnlocal(3,3) = 'NR1'
      tnlocal(3,4) = 'CPH2'
c     --- CT2 - CPH1 - NR2 - CPH2
      tnlocal(4,1) = 'CT2'
      tnlocal(4,2) = 'CPH1'
      tnlocal(4,3) = 'NR2'
      tnlocal(4,4) = 'CPH2'
c     --- CT2 - CY - CPT - CA
      tnlocal(5,1) = 'CT2'
      tnlocal(5,2) = 'CY'
      tnlocal(5,3) = 'CPT'
      tnlocal(5,4) = 'CA'
c     --- CT2 - CY - CPT - CPT
      tnlocal(6,1) = 'CT2'
      tnlocal(6,2) = 'CY'
      tnlocal(6,3) = 'CPT'
      tnlocal(6,4) = 'CPT'
c     --- CT3 - CPH1 - NR1 - CPH2
      tnlocal(7,1) = 'CT3'
      tnlocal(7,2) = 'CPH1'
      tnlocal(7,3) = 'NR1'
      tnlocal(7,4) = 'CPH2'
c     --- CY - CPT - CA - CA
      tnlocal(8,1) = 'CY'
      tnlocal(8,2) = 'CPT'
      tnlocal(8,3) = 'CA'
      tnlocal(8,4) = 'CA'
c     --- H - NR3 - CPH1 - CT2
      tnlocal(9,1) = 'H'
      tnlocal(9,2) = 'NR3'
      tnlocal(9,3) = 'CPH1'
      tnlocal(9,4) = 'CT2'
c     --- H - NR3 - CPH1 - CT3
      tnlocal(10,1) = 'H'
      tnlocal(10,2) = 'NR3'
      tnlocal(10,3) = 'CPH1'
      tnlocal(10,4) = 'CT3'
c     --- HA - CA - CPT - CPT
      tnlocal(11,1) = 'HA'
      tnlocal(11,2) = 'CA'
      tnlocal(11,3) = 'CPT'
      tnlocal(11,4) = 'CPT'
c     --- HA - CA - NY - CPT
      tnlocal(12,1) = 'HA'
      tnlocal(12,2) = 'CA'
      tnlocal(12,3) = 'NY'
      tnlocal(12,4) = 'CPT'
c     --- HA - CY - CPT - CA
      tnlocal(13,1) = 'HA'
      tnlocal(13,2) = 'CY'
      tnlocal(13,3) = 'CPT'
      tnlocal(13,4) = 'CA'
c     --- HA - CY - CPT - CPT
      tnlocal(14,1) = 'HA'
      tnlocal(14,2) = 'CY'
      tnlocal(14,3) = 'CPT'
      tnlocal(14,4) = 'CPT'
c     --- HP - CA - CA - CPT
      tnlocal(15,1) = 'HP'
      tnlocal(15,2) = 'CA'
      tnlocal(15,3) = 'CA'
      tnlocal(15,4) = 'CPT'
c     --- HP - CA - CPT - CPT
      tnlocal(16,1) = 'HP'
      tnlocal(16,2) = 'CA'
      tnlocal(16,3) = 'CPT'
      tnlocal(16,4) = 'CPT'
c     --- HP - CA - CPT - CY
      tnlocal(17,1) = 'HP'
      tnlocal(17,2) = 'CA'
      tnlocal(17,3) = 'CPT'
      tnlocal(17,4) = 'CY'
c     --- HR1 - CPH1 - NR3 - H
      tnlocal(18,1) = 'HR1'
      tnlocal(18,2) = 'CPH1'
      tnlocal(18,3) = 'NR3'
      tnlocal(18,4) = 'H'
c     --- HR1 - CPH2 - NR1 - CPH1
      tnlocal(19,1) = 'HR1'
      tnlocal(19,2) = 'CPH2'
      tnlocal(19,3) = 'NR1'
      tnlocal(19,4) = 'CPH1'
c     --- HR1 - CPH2 - NR2 - CPH1
      tnlocal(20,1) = 'HR1'
      tnlocal(20,2) = 'CPH2'
      tnlocal(20,3) = 'NR2'
      tnlocal(20,4) = 'CPH1'
c     --- HR2 - CPH2 - NR3 - CPH1
      tnlocal(21,1) = 'HR2'
      tnlocal(21,2) = 'CPH2'
      tnlocal(21,3) = 'NR3'
      tnlocal(21,4) = 'CPH1'
c     --- HR3 - CPH1 - NR1 - CPH2
      tnlocal(22,1) = 'HR3'
      tnlocal(22,2) = 'CPH1'
      tnlocal(22,3) = 'NR1'
      tnlocal(22,4) = 'CPH2'
c     --- HR3 - CPH1 - NR2 - CPH2
      tnlocal(23,1) = 'HR3'
      tnlocal(23,2) = 'CPH1'
      tnlocal(23,3) = 'NR2'
      tnlocal(23,4) = 'CPH2'
c     --- NR1 - CPH1 - CPH1 - CT2
      tnlocal(24,1) = 'NR1'
      tnlocal(24,2) = 'CPH1'
      tnlocal(24,3) = 'CPH1'
      tnlocal(24,4) = 'CT2'
c     --- NR1 - CPH1 - CPH1 - CT3
      tnlocal(25,1) = 'NR1'
      tnlocal(25,2) = 'CPH1'
      tnlocal(25,3) = 'CPH1'
      tnlocal(25,4) = 'CT3'
c     --- NR1 - CPH1 - CPH1 - HR3
      tnlocal(26,1) = 'NR1'
      tnlocal(26,2) = 'CPH1'
      tnlocal(26,3) = 'CPH1'
      tnlocal(26,4) = 'HR3'
c     --- NR2 - CPH1 - CPH1 - CT2
      tnlocal(27,1) = 'NR2'
      tnlocal(27,2) = 'CPH1'
      tnlocal(27,3) = 'CPH1'
      tnlocal(27,4) = 'CT2'
c     --- NR2 - CPH1 - CPH1 - CT3
      tnlocal(28,1) = 'NR2'
      tnlocal(28,2) = 'CPH1'
      tnlocal(28,3) = 'CPH1'
      tnlocal(28,4) = 'CT3'
c     --- NR2 - CPH1 - CPH1 - HR3
      tnlocal(29,1) = 'NR2'
      tnlocal(29,2) = 'CPH1'
      tnlocal(29,3) = 'CPH1'
      tnlocal(29,4) = 'HR3'
c     --- NY - CPT - CA - HP
      tnlocal(30,1) = 'NY'
      tnlocal(30,2) = 'CPT'
      tnlocal(30,3) = 'CA'
      tnlocal(30,4) = 'HP'
c     --- X - CPB - C - X
      tnlocal(31,1) = 'wild'
      tnlocal(31,2) = 'CPB'
      tnlocal(31,3) = 'C'
      tnlocal(31,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

C     --- Null Torsion: zero energy
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,0)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 71
c     --- CD - CP1 - N - C
      tnlocal(1,1) = 'CD'
      tnlocal(1,2) = 'CP1'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'C'
c     --- CP2 - CP3 - N - C
      tnlocal(2,1) = 'CP2'
      tnlocal(2,2) = 'CP3'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'C'
c     --- H - NH1 - CT1 - C
      tnlocal(3,1) = 'H'
      tnlocal(3,2) = 'NH1'
      tnlocal(3,3) = 'CT1'
      tnlocal(3,4) = 'C'
c     --- H - NH1 - CT1 - CC
      tnlocal(4,1) = 'H'
      tnlocal(4,2) = 'NH1'
      tnlocal(4,3) = 'CT1'
      tnlocal(4,4) = 'CC'
c     --- H - NH1 - CT1 - CD
      tnlocal(5,1) = 'H'
      tnlocal(5,2) = 'NH1'
      tnlocal(5,3) = 'CT1'
      tnlocal(5,4) = 'CD'
c     --- H - NH1 - CT1 - CT1
      tnlocal(6,1) = 'H'
      tnlocal(6,2) = 'NH1'
      tnlocal(6,3) = 'CT1'
      tnlocal(6,4) = 'CT1'
c     --- H - NH1 - CT1 - CT2
      tnlocal(7,1) = 'H'
      tnlocal(7,2) = 'NH1'
      tnlocal(7,3) = 'CT1'
      tnlocal(7,4) = 'CT2'
c     --- H - NH1 - CT1 - CT3
      tnlocal(8,1) = 'H'
      tnlocal(8,2) = 'NH1'
      tnlocal(8,3) = 'CT1'
      tnlocal(8,4) = 'CT3'
c     --- H - NH1 - CT2 - C
      tnlocal(9,1) = 'H'
      tnlocal(9,2) = 'NH1'
      tnlocal(9,3) = 'CT2'
      tnlocal(9,4) = 'C'
c     --- H - NH1 - CT2 - CC
      tnlocal(10,1) = 'H'
      tnlocal(10,2) = 'NH1'
      tnlocal(10,3) = 'CT2'
      tnlocal(10,4) = 'CC'
c     --- H - NH1 - CT2 - CD
      tnlocal(11,1) = 'H'
      tnlocal(11,2) = 'NH1'
      tnlocal(11,3) = 'CT2'
      tnlocal(11,4) = 'CD'
c     --- H - NH1 - CT2 - CT2
      tnlocal(12,1) = 'H'
      tnlocal(12,2) = 'NH1'
      tnlocal(12,3) = 'CT2'
      tnlocal(12,4) = 'CT2'
c     --- H - NH1 - CT2 - CT3
      tnlocal(13,1) = 'H'
      tnlocal(13,2) = 'NH1'
      tnlocal(13,3) = 'CT2'
      tnlocal(13,4) = 'CT3'
c     --- HA - CP3 - N - C
      tnlocal(14,1) = 'HA'
      tnlocal(14,2) = 'CP3'
      tnlocal(14,3) = 'N'
      tnlocal(14,4) = 'C'
c     --- HA - CT2 - CPH1 - CPH1
      tnlocal(15,1) = 'HA'
      tnlocal(15,2) = 'CT2'
      tnlocal(15,3) = 'CPH1'
      tnlocal(15,4) = 'CPH1'
c     --- HA - CT2 - NH1 - C
      tnlocal(16,1) = 'HA'
      tnlocal(16,2) = 'CT2'
      tnlocal(16,3) = 'NH1'
      tnlocal(16,4) = 'C'
c     --- HA - CT2 - NH1 - H
      tnlocal(17,1) = 'HA'
      tnlocal(17,2) = 'CT2'
      tnlocal(17,3) = 'NH1'
      tnlocal(17,4) = 'H'
c     --- HA - CT3 - CPH1 - CPH1
      tnlocal(18,1) = 'HA'
      tnlocal(18,2) = 'CT3'
      tnlocal(18,3) = 'CPH1'
      tnlocal(18,4) = 'CPH1'
c     --- HA - CT3 - NH1 - C
      tnlocal(19,1) = 'HA'
      tnlocal(19,2) = 'CT3'
      tnlocal(19,3) = 'NH1'
      tnlocal(19,4) = 'C'
c     --- HA - CT3 - NH1 - H
      tnlocal(20,1) = 'HA'
      tnlocal(20,2) = 'CT3'
      tnlocal(20,3) = 'NH1'
      tnlocal(20,4) = 'H'
c     --- HB - CT1 - NH1 - C
      tnlocal(21,1) = 'HB'
      tnlocal(21,2) = 'CT1'
      tnlocal(21,3) = 'NH1'
      tnlocal(21,4) = 'C'
c     --- HB - CT1 - NH1 - H
      tnlocal(22,1) = 'HB'
      tnlocal(22,2) = 'CT1'
      tnlocal(22,3) = 'NH1'
      tnlocal(22,4) = 'H'
c     --- HB - CT2 - NH1 - C
      tnlocal(23,1) = 'HB'
      tnlocal(23,2) = 'CT2'
      tnlocal(23,3) = 'NH1'
      tnlocal(23,4) = 'C'
c     --- HB - CT2 - NH1 - H
      tnlocal(24,1) = 'HB'
      tnlocal(24,2) = 'CT2'
      tnlocal(24,3) = 'NH1'
      tnlocal(24,4) = 'H'
c     --- HB - CT3 - NH1 - C
      tnlocal(25,1) = 'HB'
      tnlocal(25,2) = 'CT3'
      tnlocal(25,3) = 'NH1'
      tnlocal(25,4) = 'C'
c     --- HB - CT3 - NH1 - H
      tnlocal(26,1) = 'HB'
      tnlocal(26,2) = 'CT3'
      tnlocal(26,3) = 'NH1'
      tnlocal(26,4) = 'H'
c     --- HR2 - CPH2 - NR3 - H
      tnlocal(27,1) = 'HR2'
      tnlocal(27,2) = 'CPH2'
      tnlocal(27,3) = 'NR3'
      tnlocal(27,4) = 'H'
c     --- N - C - CT1 - CT1
      tnlocal(28,1) = 'N'
      tnlocal(28,2) = 'C'
      tnlocal(28,3) = 'CT1'
      tnlocal(28,4) = 'CT1'
c     --- N - C - CT1 - CT2
      tnlocal(29,1) = 'N'
      tnlocal(29,2) = 'C'
      tnlocal(29,3) = 'CT1'
      tnlocal(29,4) = 'CT2'
c     --- N - C - CT1 - CT3
      tnlocal(30,1) = 'N'
      tnlocal(30,2) = 'C'
      tnlocal(30,3) = 'CT1'
      tnlocal(30,4) = 'CT3'
c     --- N - C - CT1 - HB
      tnlocal(31,1) = 'N'
      tnlocal(31,2) = 'C'
      tnlocal(31,3) = 'CT1'
      tnlocal(31,4) = 'HB'
c     --- N - C - CT2 - HB
      tnlocal(32,1) = 'N'
      tnlocal(32,2) = 'C'
      tnlocal(32,3) = 'CT2'
      tnlocal(32,4) = 'HB'
c     --- N - C - CT3 - HA
      tnlocal(33,1) = 'N'
      tnlocal(33,2) = 'C'
      tnlocal(33,3) = 'CT3'
      tnlocal(33,4) = 'HA'
c     --- NH1 - C - CT1 - CT1
      tnlocal(34,1) = 'NH1'
      tnlocal(34,2) = 'C'
      tnlocal(34,3) = 'CT1'
      tnlocal(34,4) = 'CT1'
c     --- NH1 - C - CT1 - CT2
      tnlocal(35,1) = 'NH1'
      tnlocal(35,2) = 'C'
      tnlocal(35,3) = 'CT1'
      tnlocal(35,4) = 'CT2'
c     --- NH1 - C - CT1 - CT3
      tnlocal(36,1) = 'NH1'
      tnlocal(36,2) = 'C'
      tnlocal(36,3) = 'CT1'
      tnlocal(36,4) = 'CT3'
c     --- NH1 - C - CT1 - HB
      tnlocal(37,1) = 'NH1'
      tnlocal(37,2) = 'C'
      tnlocal(37,3) = 'CT1'
      tnlocal(37,4) = 'HB'
c     --- NH1 - C - CT2 - CT2
      tnlocal(38,1) = 'NH1'
      tnlocal(38,2) = 'C'
      tnlocal(38,3) = 'CT2'
      tnlocal(38,4) = 'CT2'
c     --- NH1 - C - CT2 - HA
      tnlocal(39,1) = 'NH1'
      tnlocal(39,2) = 'C'
      tnlocal(39,3) = 'CT2'
      tnlocal(39,4) = 'HA'
c     --- NH1 - C - CT2 - HB
      tnlocal(40,1) = 'NH1'
      tnlocal(40,2) = 'C'
      tnlocal(40,3) = 'CT2'
      tnlocal(40,4) = 'HB'
c     --- NH1 - C - CT3 - HA
      tnlocal(41,1) = 'NH1'
      tnlocal(41,2) = 'C'
      tnlocal(41,3) = 'CT3'
      tnlocal(41,4) = 'HA'
c     --- NH2 - CC - CT2 - HA
      tnlocal(42,1) = 'NH2'
      tnlocal(42,2) = 'CC'
      tnlocal(42,3) = 'CT2'
      tnlocal(42,4) = 'HA'
c     --- O - C - CT1 - HB
      tnlocal(43,1) = 'O'
      tnlocal(43,2) = 'C'
      tnlocal(43,3) = 'CT1'
      tnlocal(43,4) = 'HB'
c     --- O - C - CT1 - NH1
      tnlocal(44,1) = 'O'
      tnlocal(44,2) = 'C'
      tnlocal(44,3) = 'CT1'
      tnlocal(44,4) = 'NH1'
c     --- O - C - CT1 - NH3
      tnlocal(45,1) = 'O'
      tnlocal(45,2) = 'C'
      tnlocal(45,3) = 'CT1'
      tnlocal(45,4) = 'NH3'
c     --- O - C - CT2 - HA
      tnlocal(46,1) = 'O'
      tnlocal(46,2) = 'C'
      tnlocal(46,3) = 'CT2'
      tnlocal(46,4) = 'HA'
c     --- O - C - CT2 - HB
      tnlocal(47,1) = 'O'
      tnlocal(47,2) = 'C'
      tnlocal(47,3) = 'CT2'
      tnlocal(47,4) = 'HB'
c     --- O - C - CT2 - NH1
      tnlocal(48,1) = 'O'
      tnlocal(48,2) = 'C'
      tnlocal(48,3) = 'CT2'
      tnlocal(48,4) = 'NH1'
c     --- O - C - CT2 - NH3
      tnlocal(49,1) = 'O'
      tnlocal(49,2) = 'C'
      tnlocal(49,3) = 'CT2'
      tnlocal(49,4) = 'NH3'
c     --- O - C - CT3 - HA
      tnlocal(50,1) = 'O'
      tnlocal(50,2) = 'C'
      tnlocal(50,3) = 'CT3'
      tnlocal(50,4) = 'HA'
c     --- O - CC - CT2 - HA
      tnlocal(51,1) = 'O'
      tnlocal(51,2) = 'CC'
      tnlocal(51,3) = 'CT2'
      tnlocal(51,4) = 'HA'
c     --- X - CP1 - C - X
      tnlocal(52,1) = 'wild'
      tnlocal(52,2) = 'CP1'
      tnlocal(52,3) = 'C'
      tnlocal(52,4) = 'wild'
c     --- X - CP1 - CC - X
      tnlocal(53,1) = 'wild'
      tnlocal(53,2) = 'CP1'
      tnlocal(53,3) = 'CC'
      tnlocal(53,4) = 'wild'
c     --- X - CP1 - CD - X
      tnlocal(54,1) = 'wild'
      tnlocal(54,2) = 'CP1'
      tnlocal(54,3) = 'CD'
      tnlocal(54,4) = 'wild'
c     --- X - CPA - CPB - X
      tnlocal(55,1) = 'wild'
      tnlocal(55,2) = 'CPA'
      tnlocal(55,3) = 'CPB'
      tnlocal(55,4) = 'wild'
c     --- X - CPA - CPM - X
      tnlocal(56,1) = 'wild'
      tnlocal(56,2) = 'CPA'
      tnlocal(56,3) = 'CPM'
      tnlocal(56,4) = 'wild'
c     --- X - CPB - CPB - X
      tnlocal(57,1) = 'wild'
      tnlocal(57,2) = 'CPB'
      tnlocal(57,3) = 'CPB'
      tnlocal(57,4) = 'wild'
c     --- X - CPB - CT2 - X
      tnlocal(58,1) = 'wild'
      tnlocal(58,2) = 'CPB'
      tnlocal(58,3) = 'CT2'
      tnlocal(58,4) = 'wild'
c     --- X - CPB - CT3 - X
      tnlocal(59,1) = 'wild'
      tnlocal(59,2) = 'CPB'
      tnlocal(59,3) = 'CT3'
      tnlocal(59,4) = 'wild'
c     --- X - CPT - CPT - X
      tnlocal(60,1) = 'wild'
      tnlocal(60,2) = 'CPT'
      tnlocal(60,3) = 'CPT'
      tnlocal(60,4) = 'wild'
c     --- X - CT1 - CD - X
      tnlocal(61,1) = 'wild'
      tnlocal(61,2) = 'CT1'
      tnlocal(61,3) = 'CD'
      tnlocal(61,4) = 'wild'
c     --- X - CT2 - CA - X
      tnlocal(62,1) = 'wild'
      tnlocal(62,2) = 'CT2'
      tnlocal(62,3) = 'CA'
      tnlocal(62,4) = 'wild'
c     --- X - CT2 - CD - X
      tnlocal(63,1) = 'wild'
      tnlocal(63,2) = 'CT2'
      tnlocal(63,3) = 'CD'
      tnlocal(63,4) = 'wild'
c     --- X - CT2 - NC2 - X
      tnlocal(64,1) = 'wild'
      tnlocal(64,2) = 'CT2'
      tnlocal(64,3) = 'NC2'
      tnlocal(64,4) = 'wild'
c     --- X - CT3 - CA - X
      tnlocal(65,1) = 'wild'
      tnlocal(65,2) = 'CT3'
      tnlocal(65,3) = 'CA'
      tnlocal(65,4) = 'wild'
c     --- X - CT3 - CD - X
      tnlocal(66,1) = 'wild'
      tnlocal(66,2) = 'CT3'
      tnlocal(66,3) = 'CD'
      tnlocal(66,4) = 'wild'
c     --- X - CT3 - NC2 - X
      tnlocal(67,1) = 'wild'
      tnlocal(67,2) = 'CT3'
      tnlocal(67,3) = 'NC2'
      tnlocal(67,4) = 'wild'
c     --- X - FE - NPH - X
      tnlocal(68,1) = 'wild'
      tnlocal(68,2) = 'FE'
      tnlocal(68,3) = 'NPH'
      tnlocal(68,4) = 'wild'
c     --- X - FE - OM - X
      tnlocal(69,1) = 'wild'
      tnlocal(69,2) = 'FE'
      tnlocal(69,3) = 'OM'
      tnlocal(69,4) = 'wild'
c     --- X - NPH - CPA - X
      tnlocal(70,1) = 'wild'
      tnlocal(70,2) = 'NPH'
      tnlocal(70,3) = 'CPA'
      tnlocal(70,4) = 'wild'
c     --- HA - C - NH1 - CT1 (formamide cap for gramicidin)
c     --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c         `High-resolution polypeptide structure in a lamellar 
c          phase lipid environment from solid state NMR derived 
c          orientational constraints'
c         Structure, 5, 1655-1669 (1997).
c     --- T.B. Woolf personal communication 08-16-2001 (M.G. Martin)
      tnlocal(71,1) = 'HA'
      tnlocal(71,2) = 'C'
      tnlocal(71,3) = 'NH1'
      tnlocal(71,4) = 'CT1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (2.750d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.300d0)*(CON_KCALTOK)
      tclocal(5) = 4.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 10
c     --- CP1 - C - N - CP1
      tnlocal(1,1) = 'CP1'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CP1'
c     --- CP3 - N - C - CP1
      tnlocal(2,1) = 'CP3'
      tnlocal(2,2) = 'N'
      tnlocal(2,3) = 'C'
      tnlocal(2,4) = 'CP1'
c     --- CT1 - C - N - CP1
      tnlocal(3,1) = 'CT1'
      tnlocal(3,2) = 'C'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'CP1'
c     --- CT1 - C - N - CP3
      tnlocal(4,1) = 'CT1'
      tnlocal(4,2) = 'C'
      tnlocal(4,3) = 'N'
      tnlocal(4,4) = 'CP3'
c     --- CT2 - C - N - CP1
      tnlocal(5,1) = 'CT2'
      tnlocal(5,2) = 'C'
      tnlocal(5,3) = 'N'
      tnlocal(5,4) = 'CP1'
c     --- CT2 - C - N - CP3
      tnlocal(6,1) = 'CT2'
      tnlocal(6,2) = 'C'
      tnlocal(6,3) = 'N'
      tnlocal(6,4) = 'CP3'
c     --- CT3 - C - N - CP1
      tnlocal(7,1) = 'CT3'
      tnlocal(7,2) = 'C'
      tnlocal(7,3) = 'N'
      tnlocal(7,4) = 'CP1'
c     --- CT3 - C - N - CP3
      tnlocal(8,1) = 'CT3'
      tnlocal(8,2) = 'C'
      tnlocal(8,3) = 'N'
      tnlocal(8,4) = 'CP3'
c     --- O - C - N - CP1
      tnlocal(9,1) = 'O'
      tnlocal(9,2) = 'C'
      tnlocal(9,3) = 'N'
      tnlocal(9,4) = 'CP1'
c     --- O - C - N - CP3
      tnlocal(10,1) = 'O'
      tnlocal(10,2) = 'C'
      tnlocal(10,3) = 'N'
      tnlocal(10,4) = 'CP3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 8
c     --- CP2 - CP3 - N - CP1
      tnlocal(1,1) = 'CP2'
      tnlocal(1,2) = 'CP3'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CP1'
c     --- CP3 - N - CP1 - C
      tnlocal(2,1) = 'CP3'
      tnlocal(2,2) = 'N'
      tnlocal(2,3) = 'CP1'
      tnlocal(2,4) = 'C'
c     --- CP3 - N - CP1 - CC
      tnlocal(3,1) = 'CP3'
      tnlocal(3,2) = 'N'
      tnlocal(3,3) = 'CP1'
      tnlocal(3,4) = 'CC'
c     --- CP3 - N - CP1 - CP2
      tnlocal(4,1) = 'CP3'
      tnlocal(4,2) = 'N'
      tnlocal(4,3) = 'CP1'
      tnlocal(4,4) = 'CP2'
c     --- HA - CP3 - N - CP1
      tnlocal(5,1) = 'HA'
      tnlocal(5,2) = 'CP3'
      tnlocal(5,3) = 'N'
      tnlocal(5,4) = 'CP1'
c     --- HB - CP1 - N - CP3
      tnlocal(6,1) = 'HB'
      tnlocal(6,2) = 'CP1'
      tnlocal(6,3) = 'N'
      tnlocal(6,4) = 'CP3'
c     --- X - CT1 - NH3 - X
      tnlocal(7,1) = 'wild'
      tnlocal(7,2) = 'CT1'
      tnlocal(7,3) = 'NH3'
      tnlocal(7,4) = 'wild'
c     --- X - CT2 - NH3 - X
      tnlocal(8,1) = 'wild'
      tnlocal(8,2) = 'CT2'
      tnlocal(8,3) = 'NH3'
      tnlocal(8,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.080d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 14
c     --- CP2 - CP3 - NP - CP1
      tnlocal(1,1) = 'CP2'
      tnlocal(1,2) = 'CP3'
      tnlocal(1,3) = 'NP'
      tnlocal(1,4) = 'CP1'
c     --- CP3 - NP - CP1 - C
      tnlocal(2,1) = 'CP3'
      tnlocal(2,2) = 'NP'
      tnlocal(2,3) = 'CP1'
      tnlocal(2,4) = 'C'
c     --- CP3 - NP - CP1 - CC
      tnlocal(3,1) = 'CP3'
      tnlocal(3,2) = 'NP'
      tnlocal(3,3) = 'CP1'
      tnlocal(3,4) = 'CC'
c     --- CP3 - NP - CP1 - CD
      tnlocal(4,1) = 'CP3'
      tnlocal(4,2) = 'NP'
      tnlocal(4,3) = 'CP1'
      tnlocal(4,4) = 'CD'
c     --- CP3 - NP - CP1 - CP2
      tnlocal(5,1) = 'CP3'
      tnlocal(5,2) = 'NP'
      tnlocal(5,3) = 'CP1'
      tnlocal(5,4) = 'CP2'
c     --- HA - CP3 - NP - CP1
      tnlocal(6,1) = 'HA'
      tnlocal(6,2) = 'CP3'
      tnlocal(6,3) = 'NP'
      tnlocal(6,4) = 'CP1'
c     --- HB - CP1 - NP - CP3
      tnlocal(7,1) = 'HB'
      tnlocal(7,2) = 'CP1'
      tnlocal(7,3) = 'NP'
      tnlocal(7,4) = 'CP3'
c     --- HC - NP - CP1 - C
      tnlocal(8,1) = 'HC'
      tnlocal(8,2) = 'NP'
      tnlocal(8,3) = 'CP1'
      tnlocal(8,4) = 'C'
c     --- HC - NP - CP1 - CC
      tnlocal(9,1) = 'HC'
      tnlocal(9,2) = 'NP'
      tnlocal(9,3) = 'CP1'
      tnlocal(9,4) = 'CC'
c     --- HC - NP - CP1 - CD
      tnlocal(10,1) = 'HC'
      tnlocal(10,2) = 'NP'
      tnlocal(10,3) = 'CP1'
      tnlocal(10,4) = 'CD'
c     --- HC - NP - CP1 - CP2
      tnlocal(11,1) = 'HC'
      tnlocal(11,2) = 'NP'
      tnlocal(11,3) = 'CP1'
      tnlocal(11,4) = 'CP2'
c     --- HC - NP - CP1 - HB
      tnlocal(12,1) = 'HC'
      tnlocal(12,2) = 'NP'
      tnlocal(12,3) = 'CP1'
      tnlocal(12,4) = 'HB'
c     --- HC - NP - CP3 - CP2
      tnlocal(13,1) = 'HC'
      tnlocal(13,2) = 'NP'
      tnlocal(13,3) = 'CP3'
      tnlocal(13,4) = 'CP2'
c     --- HC - NP - CP3 - HA
      tnlocal(14,1) = 'HC'
      tnlocal(14,2) = 'NP'
      tnlocal(14,3) = 'CP3'
      tnlocal(14,4) = 'HA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (14.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- CPH2 - NR1 - CPH1 - CPH1
      tnlocal(1,1) = 'CPH2'
      tnlocal(1,2) = 'NR1'
      tnlocal(1,3) = 'CPH1'
      tnlocal(1,4) = 'CPH1'
c     --- CPH2 - NR2 - CPH1 - CPH1
      tnlocal(2,1) = 'CPH2'
      tnlocal(2,2) = 'NR2'
      tnlocal(2,3) = 'CPH1'
      tnlocal(2,4) = 'CPH1'
c     --- NR1 - CPH2 - NR2 - CPH1
      tnlocal(3,1) = 'NR1'
      tnlocal(3,2) = 'CPH2'
      tnlocal(3,3) = 'NR2'
      tnlocal(3,4) = 'CPH1'
c     --- NR2 - CPH1 - CPH1 - NR1
      tnlocal(4,1) = 'NR2'
      tnlocal(4,2) = 'CPH1'
      tnlocal(4,3) = 'CPH1'
      tnlocal(4,4) = 'NR1'
c     --- NR2 - CPH2 - NR1 - CPH1
      tnlocal(5,1) = 'NR2'
      tnlocal(5,2) = 'CPH2'
      tnlocal(5,3) = 'NR1'
      tnlocal(5,4) = 'CPH1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (12.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- CPH2 - NR3 - CPH1 - CPH1
      tnlocal(1,1) = 'CPH2'
      tnlocal(1,2) = 'NR3'
      tnlocal(1,3) = 'CPH1'
      tnlocal(1,4) = 'CPH1'
c     --- NR3 - CPH1 - CPH1 - NR3
      tnlocal(2,1) = 'NR3'
      tnlocal(2,2) = 'CPH1'
      tnlocal(2,3) = 'CPH1'
      tnlocal(2,4) = 'NR3'
c     --- NR3 - CPH2 - NR3 - CPH1
      tnlocal(3,1) = 'NR3'
      tnlocal(3,2) = 'CPH2'
      tnlocal(3,3) = 'NR3'
      tnlocal(3,4) = 'CPH1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- CPT - CPT - CY - CA
      tnlocal(1,1) = 'CPT'
      tnlocal(1,2) = 'CPT'
      tnlocal(1,3) = 'CY'
      tnlocal(1,4) = 'CA'
c     --- HA - CA - CPT - CY
      tnlocal(2,1) = 'HA'
      tnlocal(2,2) = 'CA'
      tnlocal(2,3) = 'CPT'
      tnlocal(2,4) = 'CY'
c     --- NY - CA - CY - CPT
      tnlocal(3,1) = 'NY'
      tnlocal(3,2) = 'CA'
      tnlocal(3,3) = 'CY'
      tnlocal(3,4) = 'CPT'
c     --- NY - CPT - CA - HA
      tnlocal(4,1) = 'NY'
      tnlocal(4,2) = 'CPT'
      tnlocal(4,3) = 'CA'
      tnlocal(4,4) = 'HA'
c     --- X - C - C - X
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'C'
      tnlocal(5,3) = 'C'
      tnlocal(5,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- CPT - CPT - NY - CA
      tnlocal(1,1) = 'CPT'
      tnlocal(1,2) = 'CPT'
      tnlocal(1,3) = 'NY'
      tnlocal(1,4) = 'CA'
c     --- CY - CA - NY - CPT
      tnlocal(2,1) = 'CY'
      tnlocal(2,2) = 'CA'
      tnlocal(2,3) = 'NY'
      tnlocal(2,4) = 'CPT'
c     --- NY - CPT - CPT - CY
      tnlocal(3,1) = 'NY'
      tnlocal(3,2) = 'CPT'
      tnlocal(3,3) = 'CPT'
      tnlocal(3,4) = 'CY'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (1.600d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (2.500d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     --- CT1 - C - NH1 - CT1
      ii = ii + 1
      tnlocal(ii,1) = 'CT1'
      tnlocal(ii,2) = 'C'
      tnlocal(ii,3) = 'NH1'
      tnlocal(ii,4) = 'CT1'
c     --- CT1 - NH1 - C - CP1
      ii = ii + 1
      tnlocal(ii,1) = 'CT1'
      tnlocal(ii,2) = 'NH1'
      tnlocal(ii,3) = 'C'
      tnlocal(ii,4) = 'CP1'
c     --- CT2 - C - NH1 - CT1
      ii = ii + 1
      tnlocal(ii,1) = 'CT2'
      tnlocal(ii,2) = 'C'
      tnlocal(ii,3) = 'NH1'
      tnlocal(ii,4) = 'CT1'
c     --- CT2 - C - NH1 - CT2
      ii = ii + 1
      tnlocal(ii,1) = 'CT2'
      tnlocal(ii,2) = 'C'
      tnlocal(ii,3) = 'NH1'
      tnlocal(ii,4) = 'CT2'
c     --- CT2 - C - NH1 - CT3
      ii = ii + 1
      tnlocal(ii,1) = 'CT2'
      tnlocal(ii,2) = 'C'
      tnlocal(ii,3) = 'NH1'
      tnlocal(ii,4) = 'CT3'
c     --- CT2 - NH1 - C - CP1
      ii = ii + 1
      tnlocal(ii,1) = 'CT2'
      tnlocal(ii,2) = 'NH1'
      tnlocal(ii,3) = 'C'
      tnlocal(ii,4) = 'CP1'
c     --- CT2 - NH1 - C - CT1
      ii = ii + 1
      tnlocal(ii,1) = 'CT2'
      tnlocal(ii,2) = 'NH1'
      tnlocal(ii,3) = 'C'
      tnlocal(ii,4) = 'CT1'
c     --- CT3 - C - NH1 - CT1
      ii = ii + 1
      tnlocal(ii,1) = 'CT3'
      tnlocal(ii,2) = 'C'
      tnlocal(ii,3) = 'NH1'
      tnlocal(ii,4) = 'CT1'
c     --- CT3 - C - NH1 - CT2
      ii = ii + 1
      tnlocal(ii,1) = 'CT3'
      tnlocal(ii,2) = 'C'
      tnlocal(ii,3) = 'NH1'
      tnlocal(ii,4) = 'CT2'
c     --- CT3 - C - NH1 - CT3
      ii = ii + 1
      tnlocal(ii,1) = 'CT3'
      tnlocal(ii,2) = 'C'
      tnlocal(ii,3) = 'NH1'
      tnlocal(ii,4) = 'CT3'
c     --- CT3 - NH1 - C - CP1
      ii = ii + 1
      tnlocal(ii,1) = 'CT3'
      tnlocal(ii,2) = 'NH1'
      tnlocal(ii,3) = 'C'
      tnlocal(ii,4) = 'CP1'
c     --- CT3 - NH1 - C - CT1
      ii = ii + 1
      tnlocal(ii,1) = 'CT3'
      tnlocal(ii,2) = 'NH1'
      tnlocal(ii,3) = 'C'
      tnlocal(ii,4) = 'CT1'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)


      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.800d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- CT1 - CT1 - NH1 - C
      tnlocal(1,1) = 'CT1'
      tnlocal(1,2) = 'CT1'
      tnlocal(1,3) = 'NH1'
      tnlocal(1,4) = 'C'
c     --- CT2 - CT1 - NH1 - C
      tnlocal(2,1) = 'CT2'
      tnlocal(2,2) = 'CT1'
      tnlocal(2,3) = 'NH1'
      tnlocal(2,4) = 'C'
c     --- CT2 - CT2 - NH1 - C
      tnlocal(3,1) = 'CT2'
      tnlocal(3,2) = 'CT2'
      tnlocal(3,3) = 'NH1'
      tnlocal(3,4) = 'C'
c     --- CT3 - CT1 - NH1 - C
      tnlocal(4,1) = 'CT3'
      tnlocal(4,2) = 'CT1'
      tnlocal(4,3) = 'NH1'
      tnlocal(4,4) = 'C'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.230d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 6
c     --- CT1 - CT2 - CA - CA
      tnlocal(1,1) = 'CT1'
      tnlocal(1,2) = 'CT2'
      tnlocal(1,3) = 'CA'
      tnlocal(1,4) = 'CA'
c     --- CT1 - CT2 - CY - CA
      tnlocal(2,1) = 'CT1'
      tnlocal(2,2) = 'CT2'
      tnlocal(2,3) = 'CY'
      tnlocal(2,4) = 'CA'
c     --- CT1 - CT2 - CY - CPT
      tnlocal(3,1) = 'CT1'
      tnlocal(3,2) = 'CT2'
      tnlocal(3,3) = 'CY'
      tnlocal(3,4) = 'CPT'
c     --- CT3 - CT2 - CA - CA
      tnlocal(4,1) = 'CT3'
      tnlocal(4,2) = 'CT2'
      tnlocal(4,3) = 'CA'
      tnlocal(4,4) = 'CA'
c     --- CT3 - CT2 - CY - CA
      tnlocal(5,1) = 'CT3'
      tnlocal(5,2) = 'CT2'
      tnlocal(5,3) = 'CY'
      tnlocal(5,4) = 'CA'
c     --- CT3 - CT2 - CY - CPT
      tnlocal(6,1) = 'CT3'
      tnlocal(6,2) = 'CT2'
      tnlocal(6,3) = 'CY'
      tnlocal(6,4) = 'CPT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.200d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.270d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.000d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CT1 - CT2 - CPH1 - CPH1
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.500d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 15
c     --- CT2 - CPH1 - NR3 - CPH2
      tnlocal(1,1) = 'CT2'
      tnlocal(1,2) = 'CPH1'
      tnlocal(1,3) = 'NR3'
      tnlocal(1,4) = 'CPH2'
c     --- H - NH1 - C - CP1
      tnlocal(2,1) = 'H'
      tnlocal(2,2) = 'NH1'
      tnlocal(2,3) = 'C'
      tnlocal(2,4) = 'CP1'
c     --- H - NH1 - C - CT1
      tnlocal(3,1) = 'H'
      tnlocal(3,2) = 'NH1'
      tnlocal(3,3) = 'C'
      tnlocal(3,4) = 'CT1'
c     --- H - NH1 - C - CT2
      tnlocal(4,1) = 'H'
      tnlocal(4,2) = 'NH1'
      tnlocal(4,3) = 'C'
      tnlocal(4,4) = 'CT2'
c     --- H - NH1 - C - CT3
      tnlocal(5,1) = 'H'
      tnlocal(5,2) = 'NH1'
      tnlocal(5,3) = 'C'
      tnlocal(5,4) = 'CT3'
c     --- H - NH2 - CC - CP1
      tnlocal(6,1) = 'H'
      tnlocal(6,2) = 'NH2'
      tnlocal(6,3) = 'CC'
      tnlocal(6,4) = 'CP1'
c     --- HA - CA - CA - HA
      tnlocal(7,1) = 'HA'
      tnlocal(7,2) = 'CA'
      tnlocal(7,3) = 'CA'
      tnlocal(7,4) = 'HA'
c     --- HR1 - CPH1 - NR3 - CPH2
      tnlocal(8,1) = 'HR1'
      tnlocal(8,2) = 'CPH1'
      tnlocal(8,3) = 'NR3'
      tnlocal(8,4) = 'CPH2'
c     --- NR3 - CPH1 - CPH1 - CT2
      tnlocal(9,1) = 'NR3'
      tnlocal(9,2) = 'CPH1'
      tnlocal(9,3) = 'CPH1'
      tnlocal(9,4) = 'CT2'
c     --- NR3 - CPH1 - CPH1 - CT3
      tnlocal(10,1) = 'NR3'
      tnlocal(10,2) = 'CPH1'
      tnlocal(10,3) = 'CPH1'
      tnlocal(10,4) = 'CT3'
c     --- NR3 - CPH1 - CPH1 - HR1
      tnlocal(11,1) = 'NR3'
      tnlocal(11,2) = 'CPH1'
      tnlocal(11,3) = 'CPH1'
      tnlocal(11,4) = 'HR1'
c     --- O - C - NH1 - CT1
      tnlocal(12,1) = 'O'
      tnlocal(12,2) = 'C'
      tnlocal(12,3) = 'NH1'
      tnlocal(12,4) = 'CT1'
c     --- O - C - NH1 - CT2
      tnlocal(13,1) = 'O'
      tnlocal(13,2) = 'C'
      tnlocal(13,3) = 'NH1'
      tnlocal(13,4) = 'CT2'
c     --- O - C - NH1 - CT3
      tnlocal(14,1) = 'O'
      tnlocal(14,2) = 'C'
      tnlocal(14,3) = 'NH1'
      tnlocal(14,4) = 'CT3'
c     --- O - C - NH1 - H
      tnlocal(15,1) = 'O'
      tnlocal(15,2) = 'C'
      tnlocal(15,3) = 'NH1'
      tnlocal(15,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.400d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 9
c     --- CT2 - CT2 - CPH1 - CPH1
      tnlocal(1,1) = 'CT2'
      tnlocal(1,2) = 'CT2'
      tnlocal(1,3) = 'CPH1'
      tnlocal(1,4) = 'CPH1'
c     --- NH1 - C - CP1 - CP2
      tnlocal(2,1) = 'NH1'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'CP1'
      tnlocal(2,4) = 'CP2'
c     --- NH1 - CT1 - C - N
      tnlocal(3,1) = 'NH1'
      tnlocal(3,2) = 'CT1'
      tnlocal(3,3) = 'C'
      tnlocal(3,4) = 'N'
c     --- NH1 - CT2 - C - N
      tnlocal(4,1) = 'NH1'
      tnlocal(4,2) = 'CT2'
      tnlocal(4,3) = 'C'
      tnlocal(4,4) = 'N'
c     --- NH3 - CT1 - C - N
      tnlocal(5,1) = 'NH3'
      tnlocal(5,2) = 'CT1'
      tnlocal(5,3) = 'C'
      tnlocal(5,4) = 'N'
c     --- NH3 - CT1 - CC - NH2
      tnlocal(6,1) = 'NH3'
      tnlocal(6,2) = 'CT1'
      tnlocal(6,3) = 'CC'
      tnlocal(6,4) = 'NH2'
c     --- NH3 - CT2 - C - N
      tnlocal(7,1) = 'NH3'
      tnlocal(7,2) = 'CT2'
      tnlocal(7,3) = 'C'
      tnlocal(7,4) = 'N'
c     --- NH3 - CT2 - C - NH1
      tnlocal(8,1) = 'NH3'
      tnlocal(8,2) = 'CT2'
      tnlocal(8,3) = 'C'
      tnlocal(8,4) = 'NH1'
c     --- NH3 - CT2 - CC - NH2
      tnlocal(9,1) = 'NH3'
      tnlocal(9,2) = 'CT2'
      tnlocal(9,3) = 'CC'
      tnlocal(9,4) = 'NH2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.150d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- CT2 - CT2 - CT2 - CT2
      tnlocal(1,1) = 'CT2'
      tnlocal(1,2) = 'CT2'
      tnlocal(1,3) = 'CT2'
      tnlocal(1,4) = 'CT2'
c     --- CT3 - CT2 - CT2 - CT2
      tnlocal(2,1) = 'CT3'
      tnlocal(2,2) = 'CT2'
      tnlocal(2,3) = 'CT2'
      tnlocal(2,4) = 'CT2'
c     --- CT3 - CT2 - CT2 - CT3
      tnlocal(3,1) = 'CT3'
      tnlocal(3,2) = 'CT2'
      tnlocal(3,3) = 'CT2'
      tnlocal(3,4) = 'CT3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (1.000d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (4.100d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.900d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- CT2 - SM - SM - CT2
      tnlocal(1,1) = 'CT2'
      tnlocal(1,2) = 'SM'
      tnlocal(1,3) = 'SM'
      tnlocal(1,4) = 'CT2'
c     --- CT3 - SM - SM - CT3
      tnlocal(2,1) = 'CT3'
      tnlocal(2,2) = 'SM'
      tnlocal(2,3) = 'SM'
      tnlocal(2,4) = 'CT3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.240d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.370d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- CT3 - CT2 - S - CT3
      tnlocal(1,1) = 'CT3'
      tnlocal(1,2) = 'CT2'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'CT3'
c     --- CT3 - S - CT2 - CT2
      tnlocal(2,1) = 'CT3'
      tnlocal(2,2) = 'S'
      tnlocal(2,3) = 'CT2'
      tnlocal(2,4) = 'CT2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (10.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- CY - CPT - CPT - CA
      tnlocal(1,1) = 'CY'
      tnlocal(1,2) = 'CPT'
      tnlocal(1,3) = 'CPT'
      tnlocal(1,4) = 'CA'
c     --- NY - CPT - CPT - CA
      tnlocal(2,1) = 'NY'
      tnlocal(2,2) = 'CPT'
      tnlocal(2,3) = 'CPT'
      tnlocal(2,4) = 'CA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.400d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 8
c     --- H - NH2 - CC - CT1
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'NH2'
      tnlocal(1,3) = 'CC'
      tnlocal(1,4) = 'CT1'
c     --- H - NH2 - CC - CT2
      tnlocal(2,1) = 'H'
      tnlocal(2,2) = 'NH2'
      tnlocal(2,3) = 'CC'
      tnlocal(2,4) = 'CT2'
c     --- H - NH2 - CC - CT3
      tnlocal(3,1) = 'H'
      tnlocal(3,2) = 'NH2'
      tnlocal(3,3) = 'CC'
      tnlocal(3,4) = 'CT3'
c     --- H - NR3 - CPH1 - CPH1
      tnlocal(4,1) = 'H'
      tnlocal(4,2) = 'NR3'
      tnlocal(4,3) = 'CPH1'
      tnlocal(4,4) = 'CPH1'
c     --- HA - CC - NH2 - H
      tnlocal(5,1) = 'HA'
      tnlocal(5,2) = 'CC'
      tnlocal(5,3) = 'NH2'
      tnlocal(5,4) = 'H'
c     --- NR3 - CPH2 - NR3 - H
      tnlocal(6,1) = 'NR3'
      tnlocal(6,2) = 'CPH2'
      tnlocal(6,3) = 'NR3'
      tnlocal(6,4) = 'H'
c     --- O - CC - NH2 - H
      tnlocal(7,1) = 'O'
      tnlocal(7,2) = 'CC'
      tnlocal(7,3) = 'NH2'
      tnlocal(7,4) = 'H'
c     --- HA - C- NH1 - H (formamide cap on gramicidin)
c     --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c         `High-resolution polypeptide structure in a lamellar 
c          phase lipid environment from solid state NMR derived 
c          orientational constraints'
c         Structure, 5, 1655-1669 (1997).
c     --- T.B. Woolf personal communication 08-16-2001 (M.G. Martin)
      tnlocal(8,1) = 'HA'
      tnlocal(8,2) = 'C'
      tnlocal(8,3) = 'NH1'
      tnlocal(8,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 11
c     --- H - NR1 - CPH1 - CPH1
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'NR1'
      tnlocal(1,3) = 'CPH1'
      tnlocal(1,4) = 'CPH1'
c     --- H - NR1 - CPH1 - CT2
      tnlocal(2,1) = 'H'
      tnlocal(2,2) = 'NR1'
      tnlocal(2,3) = 'CPH1'
      tnlocal(2,4) = 'CT2'
c     --- H - NR1 - CPH1 - CT3
      tnlocal(3,1) = 'H'
      tnlocal(3,2) = 'NR1'
      tnlocal(3,3) = 'CPH1'
      tnlocal(3,4) = 'CT3'
c     --- HA - CA  - NY - H
      tnlocal(4,1) = 'HA'
      tnlocal(4,2) = 'CA'
      tnlocal(4,3) = 'NY'
      tnlocal(4,4) = 'H'
c     --- HP - CY - CA - HP
      tnlocal(5,1) = 'HP'
      tnlocal(5,2) = 'CY'
      tnlocal(5,3) = 'CA'
      tnlocal(5,4) = 'HP'
c     --- HR1 - CPH1 - CPH1 - CT2
      tnlocal(6,1) = 'HR1'
      tnlocal(6,2) = 'CPH1'
      tnlocal(6,3) = 'CPH1'
      tnlocal(6,4) = 'CT2'
c     --- HR1 - CPH1 - CPH1 - CT3
      tnlocal(7,1) = 'HR1'
      tnlocal(7,2) = 'CPH1'
      tnlocal(7,3) = 'CPH1'
      tnlocal(7,4) = 'CT3'
c     --- HR1 - CPH1 - CPH1 - HR1
      tnlocal(8,1) = 'HR1'
      tnlocal(8,2) = 'CPH1'
      tnlocal(8,3) = 'CPH1'
      tnlocal(8,4) = 'HR1'
c     --- HR1 - CPH2 - NR1 - H
      tnlocal(9,1) = 'HR1'
      tnlocal(9,2) = 'CPH2'
      tnlocal(9,3) = 'NR1'
      tnlocal(9,4) = 'H'
c     --- HR3 - CPH1 - NR1 - H
      tnlocal(10,1) = 'HR3'
      tnlocal(10,2) = 'CPH1'
      tnlocal(10,3) = 'NR1'
      tnlocal(10,4) = 'H'
c     --- NR2 - CPH2 - NR1 - H
      tnlocal(11,1) = 'NR2'
      tnlocal(11,2) = 'CPH2'
      tnlocal(11,3) = 'NR1'
      tnlocal(11,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.800d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- H - NY - CA - CY
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'NY'
      tnlocal(1,3) = 'CA'
      tnlocal(1,4) = 'CY'
c     --- H - NY - CPT - CA
      tnlocal(2,1) = 'H'
      tnlocal(2,2) = 'NY'
      tnlocal(2,3) = 'CPT'
      tnlocal(2,4) = 'CA'
c     --- H - NY - CPT - CPT
      tnlocal(3,1) = 'H'
      tnlocal(3,2) = 'NY'
      tnlocal(3,3) = 'CPT'
      tnlocal(3,4) = 'CPT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.990d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- H - OH1 - CA - CA
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'OH1'
      tnlocal(1,3) = 'CA'
      tnlocal(1,4) = 'CA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (1.330d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.180d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.320d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- H - OH1 - CT1 - CT1
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'OH1'
      tnlocal(1,3) = 'CT1'
      tnlocal(1,4) = 'CT1'
c     --- H - OH1 - CT1 - CT3
      tnlocal(2,1) = 'H'
      tnlocal(2,2) = 'OH1'
      tnlocal(2,3) = 'CT1'
      tnlocal(2,4) = 'CT3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (1.300d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.300d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.420d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- H - OH1 - CT2 - CT1
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'OH1'
      tnlocal(1,3) = 'CT2'
      tnlocal(1,4) = 'CT1'
c     --- H - OH1 - CT2 - CT2
      tnlocal(2,1) = 'H'
      tnlocal(2,2) = 'OH1'
      tnlocal(2,3) = 'CT2'
      tnlocal(2,4) = 'CT2'
c     --- H - OH1 - CT2 - CT3
      tnlocal(3,1) = 'H'
      tnlocal(3,2) = 'OH1'
      tnlocal(3,3) = 'CT2'
      tnlocal(3,4) = 'CT3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.500d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- HA - CA - CA - CA
      tnlocal(1,1) = 'HA'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'CA'
      tnlocal(1,4) = 'CA'
c     --- HA - CA - CA - CPT
      tnlocal(2,1) = 'HA'
      tnlocal(2,2) = 'CA'
      tnlocal(2,3) = 'CA'
      tnlocal(2,4) = 'CPT'
c     --- NY - CA - CY - CT2
      tnlocal(3,1) = 'NY'
      tnlocal(3,2) = 'CA'
      tnlocal(3,3) = 'CY'
      tnlocal(3,4) = 'CT2'
c     --- NY - CA - CY - HA
      tnlocal(4,1) = 'NY'
      tnlocal(4,2) = 'CA'
      tnlocal(4,3) = 'CY'
      tnlocal(4,4) = 'HA'
c     --- NY - CA - CY - HP
      tnlocal(5,1) = 'NY'
      tnlocal(5,2) = 'CA'
      tnlocal(5,3) = 'CY'
      tnlocal(5,4) = 'HP'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.200d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- HA - CA - CY - CPT
      tnlocal(1,1) = 'HA'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'CY'
      tnlocal(1,4) = 'CPT'
c     --- HA - CA - CY - CT2
      tnlocal(2,1) = 'HA'
      tnlocal(2,2) = 'CA'
      tnlocal(2,3) = 'CY'
      tnlocal(2,4) = 'CT2'
c     --- HA - CY - CA - CPT
      tnlocal(3,1) = 'HA'
      tnlocal(3,2) = 'CY'
      tnlocal(3,3) = 'CA'
      tnlocal(3,4) = 'CPT'
c     --- HA - CY - CA - HA
      tnlocal(4,1) = 'HA'
      tnlocal(4,2) = 'CY'
      tnlocal(4,3) = 'CA'
      tnlocal(4,4) = 'HA'
c     --- HP - CA - CY - CT2
      tnlocal(5,1) = 'HP'
      tnlocal(5,2) = 'CA'
      tnlocal(5,3) = 'CY'
      tnlocal(5,4) = 'CT2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.250d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- HA - CT2 - CY - CA
      tnlocal(1,1) = 'HA'
      tnlocal(1,2) = 'CT2'
      tnlocal(1,3) = 'CY'
      tnlocal(1,4) = 'CA'
c     --- HA - CT2 - CY - CPT
      tnlocal(2,1) = 'HA'
      tnlocal(2,2) = 'CT2'
      tnlocal(2,3) = 'CY'
      tnlocal(2,4) = 'CPT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.280d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- HA - CT2 - S - CT3
      tnlocal(1,1) = 'HA'
      tnlocal(1,2) = 'CT2'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'CT3'
c     --- HA - CT3 - S - CT2
      tnlocal(2,1) = 'HA'
      tnlocal(2,2) = 'CT3'
      tnlocal(2,3) = 'S'
      tnlocal(2,4) = 'CT2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.160d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 8
c     --- HA - CT3 - CS - HA
      tnlocal(1,1) = 'HA'
      tnlocal(1,2) = 'CT3'
      tnlocal(1,3) = 'CS'
      tnlocal(1,4) = 'HA'
c     --- OC - CC - CP1 - CP2
      tnlocal(2,1) = 'OC'
      tnlocal(2,2) = 'CC'
      tnlocal(2,3) = 'CP1'
      tnlocal(2,4) = 'CP2'
c     --- OC - CC - CP1 - HB
      tnlocal(3,1) = 'OC'
      tnlocal(3,2) = 'CC'
      tnlocal(3,3) = 'CP1'
      tnlocal(3,4) = 'HB'
c     --- OC - CC - CP1 - HB
      tnlocal(4,1) = 'OC'
      tnlocal(4,2) = 'CC'
      tnlocal(4,3) = 'CP1'
      tnlocal(4,4) = 'HB'
c     --- OC - CC - CP1 - N
      tnlocal(5,1) = 'OC'
      tnlocal(5,2) = 'CC'
      tnlocal(5,3) = 'CP1'
      tnlocal(5,4) = 'N'
c     --- OC - CC - CP1 - NP 
      tnlocal(6,1) = 'OC'
      tnlocal(6,2) = 'CC'
      tnlocal(6,3) = 'CP1'
      tnlocal(6,4) = 'NP'
c     --- X - CP2 - CP2 - X
      tnlocal(7,1) = 'wild'
      tnlocal(7,2) = 'CP2'
      tnlocal(7,3) = 'CP2'
      tnlocal(7,4) = 'wild'
c     --- X - CT2 - CT3 - X
      tnlocal(8,1) = 'wild'
      tnlocal(8,2) = 'CT2'
      tnlocal(8,3) = 'CT3'
      tnlocal(8,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.200d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- HP - CA - CA - CA
      tnlocal(1,1) = 'HP'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'CA'
      tnlocal(1,4) = 'CA'
c     --- HP - CA - CA - CT2
      tnlocal(2,1) = 'HP'
      tnlocal(2,2) = 'CA'
      tnlocal(2,3) = 'CA'
      tnlocal(2,4) = 'CT2'
c     --- HP - CA - CA - CT3
      tnlocal(3,1) = 'HP'
      tnlocal(3,2) = 'CA'
      tnlocal(3,3) = 'CA'
      tnlocal(3,4) = 'CT3'
c     --- OC - CA - CA - HP
      tnlocal(4,1) = 'OC'
      tnlocal(4,2) = 'CA'
      tnlocal(4,3) = 'CA'
      tnlocal(4,4) = 'HP'
c     --- OH1 - CA - CA - HP
      tnlocal(5,1) = 'OH1'
      tnlocal(5,2) = 'CA'
      tnlocal(5,3) = 'CA'
      tnlocal(5,4) = 'HP'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.400d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HP - CA - CA - HP
      tnlocal(1,1) = 'HP'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'CA'
      tnlocal(1,4) = 'HP'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- HP - CA - CY - CPT
      tnlocal(1,1) = 'HP'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'CY'
      tnlocal(1,4) = 'CPT'
c     --- HP - CA - NY - CPT
      tnlocal(2,1) = 'HP'
      tnlocal(2,2) = 'CA'
      tnlocal(2,3) = 'NY'
      tnlocal(2,4) = 'CPT'
c     --- HR3 - CPH1 - CPH1 - CT2
      tnlocal(3,1) = 'HR3'
      tnlocal(3,2) = 'CPH1'
      tnlocal(3,3) = 'CPH1'
      tnlocal(3,4) = 'CT2'
c     --- HR3 - CPH1 - CPH1 - CT3
      tnlocal(4,1) = 'HR3'
      tnlocal(4,2) = 'CPH1'
      tnlocal(4,3) = 'CPH1'
      tnlocal(4,4) = 'CT3'
c     --- HR3 - CPH1 - CPH1 - HR3
      tnlocal(5,1) = 'HR3'
      tnlocal(5,2) = 'CPH1'
      tnlocal(5,3) = 'CPH1'
      tnlocal(5,4) = 'HR3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.400d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HP - CA - NY - H
      tnlocal(1,1) = 'HP'
      tnlocal(1,2) = 'CA'
      tnlocal(1,3) = 'NY'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.800d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- HP - CY - CPT - CA
      tnlocal(1,1) = 'HP'
      tnlocal(1,2) = 'CY'
      tnlocal(1,3) = 'CPT'
      tnlocal(1,4) = 'CA'
c     --- HP - CY - CPT - CPT
      tnlocal(2,1) = 'HP'
      tnlocal(2,2) = 'CY'
      tnlocal(2,3) = 'CPT'
      tnlocal(2,4) = 'CPT'
c     --- NY - CPT - CA - CA
      tnlocal(3,1) = 'NY'
      tnlocal(3,2) = 'CPT'
      tnlocal(3,3) = 'CA'
      tnlocal(3,4) = 'CA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.240d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.150d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.270d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- HS - S - CT2 - CT1
      tnlocal(1,1) = 'HS'
      tnlocal(1,2) = 'S'
      tnlocal(1,3) = 'CT2'
      tnlocal(1,4) = 'CT1'
c     --- HS - S - CT2 - CT3
      tnlocal(2,1) = 'HS'
      tnlocal(2,2) = 'S'
      tnlocal(2,3) = 'CT2'
      tnlocal(2,4) = 'CT3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.200d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- HS - S - CT2 - HA
      tnlocal(1,1) = 'HS'
      tnlocal(1,2) = 'S'
      tnlocal(1,3) = 'CT2'
      tnlocal(1,4) = 'HA'
c     --- HS - S - CT3 - HA
      tnlocal(2,1) = 'HS'
      tnlocal(2,2) = 'S'
      tnlocal(2,3) = 'CT3'
      tnlocal(2,4) = 'HA'
c     --- X - CT1 - CT1 - X
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CT1'
      tnlocal(3,3) = 'CT1'
      tnlocal(3,4) = 'wild'
c     --- X - CT1 - CT2 - X
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'CT1'
      tnlocal(4,3) = 'CT2'
      tnlocal(4,4) = 'wild'
c     --- X - CT1 - CT3 - X
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'CT1'
      tnlocal(5,3) = 'CT3'
      tnlocal(5,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.400d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.600d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- N - C - CP1 - HB
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CP1'
      tnlocal(1,4) = 'HB'
c     --- NH1 - C - CP1 - HB
      tnlocal(2,1) = 'NH1'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'CP1'
      tnlocal(2,4) = 'HB'
c     --- NH2 - CC - CP1 - HB
      tnlocal(3,1) = 'NH2'
      tnlocal(3,2) = 'CC'
      tnlocal(3,3) = 'CP1'
      tnlocal(3,4) = 'HB'
c     --- O - C - CP1 - CP2
      tnlocal(4,1) = 'O'
      tnlocal(4,2) = 'C'
      tnlocal(4,3) = 'CP1'
      tnlocal(4,4) = 'CP2'
c     --- O - CC - CP1 - CP2
      tnlocal(5,1) = 'O'
      tnlocal(5,2) = 'CC'
      tnlocal(5,3) = 'CP1'
      tnlocal(5,4) = 'CP2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.300d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (-0.300d0)*(CON_KCALTOK)
      tclocal(5) = 4.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- N - C - CP1 - N
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CP1'
      tnlocal(1,4) = 'N'
c     --- NH1 - C - CP1 - N
      tnlocal(2,1) = 'NH1'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'CP1'
      tnlocal(2,4) = 'N'
c     --- NH2 - CC - CP1 - N
      tnlocal(3,1) = 'NH2'
      tnlocal(3,2) = 'CC'
      tnlocal(3,3) = 'CP1'
      tnlocal(3,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.600d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- NH1 - C - CT1 - NH1
      tnlocal(1,1) = 'NH1'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CT1'
      tnlocal(1,4) = 'NH1'
c     --- NH1 - C - CT2 - NH1
      tnlocal(2,1) = 'NH1'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'CT2'
      tnlocal(2,4) = 'NH1'
c     --- NH3 - CT1 - C - NH1
      tnlocal(3,1) = 'NH3'
      tnlocal(3,2) = 'CT1'
      tnlocal(3,3) = 'C'
      tnlocal(3,4) = 'NH1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.400d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.600d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- N - C - CP1 - CP2
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CP1'
      tnlocal(1,4) = 'CP2'
c     --- NH1 - C - CP1 - CP2
      tnlocal(2,1) = 'NH1'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'CP1'
      tnlocal(2,4) = 'CP2'
c     --- NH2 - CC - CP1 - CP2
      tnlocal(3,1) = 'NH2'
      tnlocal(3,2) = 'CC'
      tnlocal(3,3) = 'CP1'
      tnlocal(3,4) = 'CP2'
c     --- O - C - CP1 - HB
      tnlocal(4,1) = 'O'
      tnlocal(4,2) = 'C'
      tnlocal(4,3) = 'CP1'
      tnlocal(4,4) = 'HB'
c     --- O - CC - CP1 - HB
      tnlocal(5,1) = 'O'
      tnlocal(5,2) = 'CC'
      tnlocal(5,3) = 'CP1'
      tnlocal(5,4) = 'HB'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.300d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- NP - CP1 - C - N
      tnlocal(1,1) = 'NP'
      tnlocal(1,2) = 'CP1'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'N'
c     --- NP - CP1 - C - NH1
      tnlocal(2,1) = 'NP'
      tnlocal(2,2) = 'CP1'
      tnlocal(2,3) = 'C'
      tnlocal(2,4) = 'NH1'
c     --- NP - CP1 - C - NH2
      tnlocal(3,1) = 'NP'
      tnlocal(3,2) = 'CP1'
      tnlocal(3,3) = 'C'
      tnlocal(3,4) = 'NH2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.190d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 15
c     --- NR1 - CPH1 - CT2 - CT1
      tnlocal(1,1) = 'NR1'
      tnlocal(1,2) = 'CPH1'
      tnlocal(1,3) = 'CT2'
      tnlocal(1,4) = 'CT1'
c     --- NR1 - CPH1 - CT2 - CT2
      tnlocal(2,1) = 'NR1'
      tnlocal(2,2) = 'CPH1'
      tnlocal(2,3) = 'CT2'
      tnlocal(2,4) = 'CT2'
c     --- NR1 - CPH1 - CT2 - CT3
      tnlocal(3,1) = 'NR1'
      tnlocal(3,2) = 'CPH1'
      tnlocal(3,3) = 'CT2'
      tnlocal(3,4) = 'CT3'
c     --- NR1 - CPH1 - CT2 - HA
      tnlocal(4,1) = 'NR1'
      tnlocal(4,2) = 'CPH1'
      tnlocal(4,3) = 'CT2'
      tnlocal(4,4) = 'HA'
c     --- NR1 - CPH1 - CT3 - HA
      tnlocal(5,1) = 'NR1'
      tnlocal(5,2) = 'CPH1'
      tnlocal(5,3) = 'CT3'
      tnlocal(5,4) = 'HA'
c     --- NR2 - CPH1 - CT2 - CT1
      tnlocal(6,1) = 'NR2'
      tnlocal(6,2) = 'CPH1'
      tnlocal(6,3) = 'CT2'
      tnlocal(6,4) = 'CT1'
c     --- NR2 - CPH1 - CT2 - CT2
      tnlocal(7,1) = 'NR2'
      tnlocal(7,2) = 'CPH1'
      tnlocal(7,3) = 'CT2'
      tnlocal(7,4) = 'CT2'
c     --- NR2 - CPH1 - CT2 - CT3
      tnlocal(8,1) = 'NR2'
      tnlocal(8,2) = 'CPH1'
      tnlocal(8,3) = 'CT2'
      tnlocal(8,4) = 'CT3'
c     --- NR2 - CPH1 - CT2 - HA
      tnlocal(9,1) = 'NR2'
      tnlocal(9,2) = 'CPH1'
      tnlocal(9,3) = 'CT2'
      tnlocal(9,4) = 'HA'
c     --- NR2 - CPH1 - CT3 - HA
      tnlocal(10,1) = 'NR2'
      tnlocal(10,2) = 'CPH1'
      tnlocal(10,3) = 'CT3'
      tnlocal(10,4) = 'HA'
c     --- NR3 - CPH1 - CT2 - CT1
      tnlocal(11,1) = 'NR3'
      tnlocal(11,2) = 'CPH1'
      tnlocal(11,3) = 'CT2'
      tnlocal(11,4) = 'CT1'
c     --- NR3 - CPH1 - CT2 - CT2
      tnlocal(12,1) = 'NR3'
      tnlocal(12,2) = 'CPH1'
      tnlocal(12,3) = 'CT2'
      tnlocal(12,4) = 'CT2'
c     --- NR3 - CPH1 - CT2 - CT3
      tnlocal(13,1) = 'NR3'
      tnlocal(13,2) = 'CPH1'
      tnlocal(13,3) = 'CT2'
      tnlocal(13,4) = 'CT3'
c     --- NR3 - CPH1 - CT2 - HA
      tnlocal(14,1) = 'NR3'
      tnlocal(14,2) = 'CPH1'
      tnlocal(14,3) = 'CT2'
      tnlocal(14,4) = 'HA'
c     --- NR3 - CPH1 - CT3 - HA
      tnlocal(15,1) = 'NR3'
      tnlocal(15,2) = 'CPH1'
      tnlocal(15,3) = 'CT3'
      tnlocal(15,4) = 'HA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (-0.300d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- O - C - CP1 - N
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CP1'
      tnlocal(1,4) = 'N'
c     --- O - CC - CP1 - N
      tnlocal(2,1) = 'O'
      tnlocal(2,2) = 'CC'
      tnlocal(2,3) = 'CP1'
      tnlocal(2,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.400d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- O - C - CT1 - CT1
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CT1'
      tnlocal(1,4) = 'CT1'
c     --- O - C - CT1 - CT2
      tnlocal(2,1) = 'O'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'CT1'
      tnlocal(2,4) = 'CT2'
c     --- O - C - CT1 - CT3
      tnlocal(3,1) = 'O'
      tnlocal(3,2) = 'C'
      tnlocal(3,3) = 'CT1'
      tnlocal(3,4) = 'CT3'
c     --- O - C - CT2 - CT2
      tnlocal(4,1) = 'O'
      tnlocal(4,2) = 'C'
      tnlocal(4,3) = 'CT2'
      tnlocal(4,4) = 'CT2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.965d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (3.850d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- OB - CD - OS - CT2
      tnlocal(1,1) = 'OB'
      tnlocal(1,2) = 'CD'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT2'
c     --- OB - CD - OS - CT3
      tnlocal(2,1) = 'OB'
      tnlocal(2,2) = 'CD'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'CT3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.200d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- OC - CC - CT1 - NH3
      tnlocal(1,1) = 'OC'
      tnlocal(1,2) = 'CC'
      tnlocal(1,3) = 'CT1'
      tnlocal(1,4) = 'NH3'
c     --- OC - CC - CT2 - NH3
      tnlocal(2,1) = 'OC'
      tnlocal(2,2) = 'CC'
      tnlocal(2,3) = 'CT2'
      tnlocal(2,4) = 'NH3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.010d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- S - CT2 - CT2 - HA
      tnlocal(1,1) = 'S'
      tnlocal(1,2) = 'CT2'
      tnlocal(1,3) = 'CT2'
      tnlocal(1,4) = 'HA'
c     --- SM - CT2 - CT2 - HA
      tnlocal(2,1) = 'SM'
      tnlocal(2,2) = 'CT2'
      tnlocal(2,3) = 'CT2'
      tnlocal(2,4) = 'HA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.310d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- SM - SM - CT2 - CT1
      tnlocal(1,1) = 'SM'
      tnlocal(1,2) = 'SM'
      tnlocal(1,3) = 'CT2'
      tnlocal(1,4) = 'CT1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.158d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- SM - SM - CT2 - HA
      tnlocal(1,1) = 'SM'
      tnlocal(1,2) = 'SM'
      tnlocal(1,3) = 'CT2'
      tnlocal(1,4) = 'HA'
c     --- SM - SM - CT3 - HA
      tnlocal(2,1) = 'SM'
      tnlocal(2,2) = 'SM'
      tnlocal(2,3) = 'CT3'
      tnlocal(2,4) = 'HA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.150d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- SS - CS - CT3 - HA
      tnlocal(1,1) = 'SS'
      tnlocal(1,2) = 'CS'
      tnlocal(1,3) = 'CT3'
      tnlocal(1,4) = 'HA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.250d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C - NC2 - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'NC2'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.050d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- X - CD - OH1 - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CD'
      tnlocal(1,3) = 'OH1'
      tnlocal(1,4) = 'wild'
c     --- X - CD - OS - X
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CD'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.140d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- X - CP1 - CP2 - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CP1'
      tnlocal(1,3) = 'CP2'
      tnlocal(1,4) = 'wild'
c     --- X - CP3 - CP2 - X
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CP3'
      tnlocal(2,3) = 'CP2'
      tnlocal(2,4) = 'wild'
c     --- X - CT1 - OH1 - X
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CT1'
      tnlocal(3,3) = 'OH1'
      tnlocal(3,4) = 'wild'
c     --- X - CT2 - OH1 - X
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'CT2'
      tnlocal(4,3) = 'OH1'
      tnlocal(4,4) = 'wild'
c     --- X - CT3 - OH1 - X
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'CT3'
      tnlocal(5,3) = 'OH1'
      tnlocal(5,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.050d0)*(CON_KCALTOK)
      tclocal(2) = 6.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- X - CT1 - CC - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT1'
      tnlocal(1,3) = 'CC'
      tnlocal(1,4) = 'wild'
c     --- X - CT2 - CC - X
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CT2'
      tnlocal(2,3) = 'CC'
      tnlocal(2,4) = 'wild'
c     --- X - CT3 - CC - X
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CT3'
      tnlocal(3,3) = 'CC'
      tnlocal(3,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (-0.100d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- X - CT1 - OS - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT1'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'wild'
c     --- X - CT2 - OS - X
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CT2'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'wild'
c     --- X - CT3 - OS - X
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CT3'
      tnlocal(3,3) = 'OS'
      tnlocal(3,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.195d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CT2 - CT2 - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT2'
      tnlocal(1,3) = 'CT2'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.155d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CT3 - CT3 - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT3'
      tnlocal(1,3) = 'CT3'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.110d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CT3 - NH2 - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT3'
      tnlocal(1,3) = 'NH2'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.090d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CT3 - NH3 - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT3'
      tnlocal(1,3) = 'NH3'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.050d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - FE - CM - X
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'FE'
      tnlocal(1,3) = 'CM'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)
      
      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Charmm22 improper harmonic torsions                            *
c     ******************************************************************
c     --- charmm units are kcal/mol/rad^2 and degrees
c     --- listed as i-j-k-l
c     --- It is not clear to me exactly what the order corresponds to.
c     --- in some cases we have l bonded to i,j,k as the central planar
c     --- atom
c     --- in other cases i is the central atom bonded to j,k,l
c     --- in either of these cases if i is the central atom then 
c     --- l is bonded only to i and not to any other atom.
c     --- Similarily if l is the central atom the i is only bonded 
c     --- to l and not to any other atom.

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (20.80d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- CPB - CPA - NPH - CPA
      ninlocal(1,1) = 'CPB'
      ninlocal(1,2) = 'CPA'
      ninlocal(1,3) = 'NPH'
      ninlocal(1,4) = 'CPA'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (90.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 3
c     --- CPB - X - X - C
      ninlocal(1,1) = 'CPB'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'C'
c     --- CT2 - X - X - CPB
      ninlocal(2,1) = 'CT2'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'CPB'
c     --- CT3 - X - X - CPB
      ninlocal(3,1) = 'CT3'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'CPB'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (20.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 3
c     --- HA - C - C - HA
      ninlocal(1,1) = 'HA'
      ninlocal(1,2) = 'C'
      ninlocal(1,3) = 'C'
      ninlocal(1,4) = 'HA'
c     --- HA - CPB - C - C
      ninlocal(2,1) = 'HA'
      ninlocal(2,2) = 'CPB'
      ninlocal(2,3) = 'C'
      ninlocal(2,4) = 'C'
c     --- NH1 - X - X - H
      ninlocal(3,1) = 'NH1'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'H'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (29.40d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- HA - CPA - CPA - CPM
      ninlocal(1,1) = 'HA'
      ninlocal(1,2) = 'CPA'
      ninlocal(1,3) = 'CPA'
      ninlocal(1,4) = 'CPM'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (20.00d0)*(CON_KCALTOK)
      iclocal(1) = (180.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- HA - HA - C - C
      ninlocal(1,1) = 'HA'
      ninlocal(1,2) = 'HA'
      ninlocal(1,3) = 'C'
      ninlocal(1,4) = 'C'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.50d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 6
c     --- HR1 - NR1 - NR2 - CPH2
      ninlocal(1,1) = 'HR1'
      ninlocal(1,2) = 'NR1'
      ninlocal(1,3) = 'NR2'
      ninlocal(1,4) = 'CPH2'
c     --- HR1 - NR2 - NR1 - CPH2
      ninlocal(2,1) = 'HR1'
      ninlocal(2,2) = 'NR2'
      ninlocal(2,3) = 'NR1'
      ninlocal(2,4) = 'CPH2'
c     --- HR3 - CPH1 - NR1 - CPH1
      ninlocal(3,1) = 'HR3'
      ninlocal(3,2) = 'CPH1'
      ninlocal(3,3) = 'NR1'
      ninlocal(3,4) = 'CPH1'
c     --- HR3 - CPH1 - NR2 - CPH1
      ninlocal(4,1) = 'HR3'
      ninlocal(4,2) = 'CPH1'
      ninlocal(4,3) = 'NR2'
      ninlocal(4,4) = 'CPH1'
c     --- HR3 - NR1 - CPH1 - CPH1
      ninlocal(5,1) = 'HR3'
      ninlocal(5,2) = 'NR1'
      ninlocal(5,3) = 'CPH1'
      ninlocal(5,4) = 'CPH1'
c     --- HR3 - NR2 - CPH1 - CPH1
      ninlocal(6,1) = 'HR3'
      ninlocal(6,2) = 'NR2'
      ninlocal(6,3) = 'CPH1'
      ninlocal(6,4) = 'CPH1'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (1.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- HR3 - CPH1 - NR3 - CPH1
      ninlocal(1,1) = 'HR3'
      ninlocal(1,2) = 'CPH1'
      ninlocal(1,3) = 'NR3'
      ninlocal(1,4) = 'CPH1'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
c     --- Null improper torsion: not clear why this is needed
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- N - C - CP1 - CP3
      ninlocal(1,1) = 'N'
      ninlocal(1,2) = 'C'
      ninlocal(1,3) = 'CP1'
      ninlocal(1,4) = 'CP3'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (40.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- NC2 - X - X - C
      ninlocal(1,1) = 'NC2'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'C'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (4.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- NH2 - X - X- H
      ninlocal(1,1) = 'NH2'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'H'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (137.40d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- NPH - CPA - CPA - FE
      ninlocal(1,1) = 'NPH'
      ninlocal(1,2) = 'CPA'
      ninlocal(1,3) = 'CPA'
      ninlocal(1,4) = 'FE'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (40.60d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- NPH - CPA - CPB - CPB
      ninlocal(1,1) = 'NPH'
      ninlocal(1,2) = 'CPA'
      ninlocal(1,3) = 'CPB'
      ninlocal(1,4) = 'CPB'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (18.30d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- NPH - CPA - CPM - CPA
      ninlocal(1,1) = 'NPH'
      ninlocal(1,2) = 'CPA'
      ninlocal(1,3) = 'CPM'
      ninlocal(1,4) = 'CPA'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (32.70d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- NPH - CPM - CPB - CPA
      ninlocal(1,1) = 'NPH'
      ninlocal(1,2) = 'CPM'
      ninlocal(1,3) = 'CPB'
      ninlocal(1,4) = 'CPA'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.45d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 2
c     --- NR1 - CPH1 - CPH2 - H
      ninlocal(1,1) = 'NR1'
      ninlocal(1,2) = 'CPH1'
      ninlocal(1,3) = 'CPH2'
      ninlocal(1,4) = 'H'
c     --- NR1 - CPH2 - CPH1 - H
      ninlocal(2,1) = 'NR1'
      ninlocal(2,2) = 'CPH2'
      ninlocal(2,3) = 'CPH1'
      ninlocal(2,4) = 'H'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (1.20d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 2
c     --- NR3 - CPH1 - CPH2 - H
      ninlocal(1,1) = 'NR3'
      ninlocal(1,2) = 'CPH1'
      ninlocal(1,3) = 'CPH2'
      ninlocal(1,4) = 'H'
c     --- NR3 - CPH2 - CPH1 - H
      ninlocal(2,1) = 'NR3'
      ninlocal(2,2) = 'CPH2'
      ninlocal(2,3) = 'CPH1'
      ninlocal(2,4) = 'H'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (100.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 2
c     --- NY - CA - CY - CPT
      ninlocal(1,1) = 'NY'
      ninlocal(1,2) = 'CA'
      ninlocal(1,3) = 'CY'
      ninlocal(1,4) = 'CPT'
c     --- OB - X - X - CD
      ninlocal(2,1) = 'OB'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'CD'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (45.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 10
c     --- O - CP1 - NH2 - CC
      ninlocal(1,1) = 'O'
      ninlocal(1,2) = 'CP1'
      ninlocal(1,3) = 'NH2'
      ninlocal(1,4) = 'CC'
c     --- O - CT1 - NH2 - CC
      ninlocal(2,1) = 'O'
      ninlocal(2,2) = 'CT1'
      ninlocal(2,3) = 'NH2'
      ninlocal(2,4) = 'CC'
c     --- O - CT2 - NH2 - CC
      ninlocal(3,1) = 'O'
      ninlocal(3,2) = 'CT2'
      ninlocal(3,3) = 'NH2'
      ninlocal(3,4) = 'CC'
c     --- O - CT3 - NH2 - CC
      ninlocal(4,1) = 'O'
      ninlocal(4,2) = 'CT3'
      ninlocal(4,3) = 'NH2'
      ninlocal(4,4) = 'CC'
c     --- O - HA - NH2 - CC
      ninlocal(5,1) = 'O'
      ninlocal(5,2) = 'HA'
      ninlocal(5,3) = 'NH2'
      ninlocal(5,4) = 'CC'
c     --- O - NH2 - CP1 - CC
      ninlocal(6,1) = 'O'
      ninlocal(6,2) = 'NH2'
      ninlocal(6,3) = 'CP1'
      ninlocal(6,4) = 'CC'
c     --- O - NH2 - CT1 - CC
      ninlocal(7,1) = 'O'
      ninlocal(7,2) = 'NH2'
      ninlocal(7,3) = 'CT1'
      ninlocal(7,4) = 'CC'
c     --- O - NH2 - CT2 - CC
      ninlocal(8,1) = 'O'
      ninlocal(8,2) = 'NH2'
      ninlocal(8,3) = 'CT2'
      ninlocal(8,4) = 'CC'
c     --- O - NH2 - CT3 - CC
      ninlocal(9,1) = 'O'
      ninlocal(9,2) = 'NH2'
      ninlocal(9,3) = 'CT3'
      ninlocal(9,4) = 'CC'
c     --- O - NH2 - HA - CC
      ninlocal(10,1) = 'O'
      ninlocal(10,2) = 'NH2'
      ninlocal(10,3) = 'HA'
      ninlocal(10,4) = 'CC'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (120.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 2
c     --- O - N - CT2 - CC
      ninlocal(1,1) = 'O'
      ninlocal(1,2) = 'N'
      ninlocal(1,3) = 'CT2'
      ninlocal(1,4) = 'CC'
c     --- O - X - X - C
      ninlocal(2,1) = 'O'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'C'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (96.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- OC - X - X - CC
      ninlocal(1,1) = 'OC'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'CC'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return
      
      return
      end
