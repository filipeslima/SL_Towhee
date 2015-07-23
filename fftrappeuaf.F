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
      subroutine twh_fftrappeuaf(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the trappe-uaf   *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 11-11-2005 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nbitypes
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
      double precision qbivlocal 
      character*(FFNAMELEN) qbinlocal 
      dimension qbinlocal(2)

      ffname = 'TraPPE-UAf'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

c     ******************************************************************
c     * TraPPE-UAf                                                     *
c     * an implementation of the TraPPE-UA force field using harmonic  *
c     * bonds taken from other work                                    *
c     ******************************************************************
c     References:
c     [TraPPE-1] 
c     M.G. Martin, J.I. Siepmann, 
c     `Transferable potentials for phase equilibria. 1. 
c     United-atom description of n-alkanes'
c     J. Phys. Chem. B 102, 2569-2577 (1998)
c     [TraPPE-2] 
c     M.G. Martin, J.I. Siepmann, 
c     `Novel configurational-bias Monte Carlo method for
c     branched molecules. Transferable potentials for 
c     phase equilibria. 2. united-atom description of
c     branched alkanes'
c     J. Phys. Chem. B 103, 4508-4517 (1999)
c     [Mundy et al 1996]
c     C.J. Mundy, S. Balasubramanian, K. Bagchi, J.I. Siepmann,
c     M.L. Klein,
c     `Equilibrium and non-equilibrium studies of fluid alkanes in bulk
c     and at interfaces', Faraday Discuss. 104, 17-36 (1996).

      type = 0

c     *** Carbon (plus bonded hydrogens) ***

c     *** CH4 molecules ***
      
c     (CH4) Methane CH4 (sp3)  [TraPPE-1]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.73d0
      nbclocal(2) = 148.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0 + 4.0d0*1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH4'
      cbnlocal = 'none'
      canlocal = 'none'
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** CH3 united-atoms ***
      
c     CH3*(sp3) Alkane methyl group (sp3) bonded to anything 
c     -------- [TraPPE-1]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75d0
      nbclocal(2) = 98.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0 + 3.0d0*1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3*(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     *** CH2 united-atoms ***

c     CH2**(sp3) Alkane [CH2] (sp3) bonded to anything [TraPPE-1]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.95d0
      nbclocal(2) = 46.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0 + 2.0d0*1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2**(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'CH2sp3'
      ctnlocal = 'CH2sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** CH united-atoms ***
      
c     CHccc(sp3) Alkane [CH] (sp3) bonded only to Carbons [TraPPE-2]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.68d0
      nbclocal(2) = 10.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0 + 1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHccc(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'CHsp3'
      ctnlocal = 'CHsp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** C atoms ***
      
c     Ccccc(sp3) Alkane [C] (sp3) bonded only to Carbons [TraPPE-2]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 6.4d0
      nbclocal(2) = 0.5d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ccccc(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'Csp3'
      ctnlocal = 'Csp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * TraPPE flexible bond lengths from Mundy 1996                   *
c     ******************************************************************

c     TraPPE fixed bond length [TraPPE-1]
c     [CHx](sp3) -- [CHy](sp3,sp2,aro) single bonds
c     --- from table 1.6 of McMurry Organic Chemistry 3rd edition
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.54d0
      vclocal(1) = (452900.0d0)*(0.5d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     Csp3 - Csp3 [TraPPE-1]
      isame = isame + 1
      vnlocal(isame,1) = 'Csp3'
      vnlocal(isame,2) = 'Csp3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * TraPPE-UA angles - simple harmonic                             *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.0d0
      bclocal(1) = (62500.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2(sp3) - CHy [TraPPE-1]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CH2sp3'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.0d0
      bclocal(1) = (62500.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH(sp3) - CHy [TraPPE-2]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CHsp3'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47d0
      bclocal(1) = (62500.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - C(sp3) - CHy [TraPPE-2]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Csp3'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * TraPPE torsions using standard OPLS cosine series              *
c     * k0 + k1*(1 + Cos(phi)) + k2*(1 - Cos(2*phi))                   *
c     * + k3*(1 + Cos(3*phi)                                           *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 355.03d0
      tclocal(2) = -68.19d0
      tclocal(3) = 791.32d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2(sp3) - CH2(sp3) - CHy [TraPPE-1]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CH2sp3'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     -- CH split OPLS torsion Siepmann et al Mol Phys 90 687-693 (1997)
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = -251.06d0
      tclocal(1) =  428.73d0
      tclocal(2) = -111.85d0
      tclocal(3) =  441.27d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2(sp3) - CH(sp3) - CHy [TraPPE 2]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CHsp3'
      tnlocal(isame,4) = 'Cwild'
c     CHx - CH(sp3) - CH(sp3) - CHy [TraPPE 2]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CHsp3'
      tnlocal(isame,3) = 'CHsp3'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     --- Mundy et al, Faraday Disc 104, 17-36 (1996)
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) =    0.00d0
      tclocal(1) =    0.00d0
      tclocal(2) =    0.00d0
      tclocal(3) =  461.29d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- torsional parameters for segment containing a quaternary 
c     --- carbon as a mid-segment 
c     CHx - CH2(sp3) - C(sp3) - CHy [TraPPE 2]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'Csp3'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * TraPPE-UA bond increments                                      *
c     * inferred from the charge examples in the various papers        *
c     ******************************************************************

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [TraPPE-1]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Csp3'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
