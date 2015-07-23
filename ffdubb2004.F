#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2014 Marcus G. Martin                       *
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
      subroutine twh_ffdubb2004(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the dubb2004     *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 12-10-2004 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nbitypes
c     --- local variables
      integer type,j,ij,isame
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

      ffname = 'Dubb2004'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Explicit')
      type = 0

c     --- references for Dubb2004
c     [prl2004]
c     D. Dubbeldam; S. Calero; T.J.H. Vlugt; R. Krishna; T.L.M Maesen;
c     E. Beerdsen; B. Smit;
c     "Force Field Parametrization through Fitting on Inflection Points
c     in Isotherms", Phys. Rev. Lett. 93, ?? (2004).
c     [jpcb2004]
c     D. Dubbeldam; S. Calero; T.J.H. Vlugt; R. Krishna; T.L.M. Maeson;
c     B. Smit;
c     "United Atom Force Field for Alkanes in Nanoporous Materials",
c     J. Phys. Chem. B 108, 12301-12313 (2004).
c     [jacs2004]
c     S. Calero; D. Dubbeldam; R. Krishna; B. Smit; T.J.H. Vlugt;
c     J.F.M. Denayer; J.A. Mertens; T.L.M. Maesen;
c     "Understanding the Role of Sodium during Adsorption: A Force Field
c     for Alkanes in the Sodium-Exchanged Faujasites",
c     J. Am. Chem. Soc. 126, 11377-11386 (2004).

c     --- we are not using the 1-4 nonbonded terms so set them to zero
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0

c     Oxygen zeolite atom 
      type = 1
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'ion'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of O with O (1-1) is not available as no value is 
c     --- listed in the paper [jpcb2004]
c     --- cross term of O with CH4 (1-2) [jpcb2004]
      j = 2
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.47d0
      nbclocal(2) = 115.00d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with CH3 (1-3) [jpcb2004]
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.48d0
      nbclocal(2) = 93.00d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with CH2 (1-4) [jpcb2004]
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.58d0
      nbclocal(2) = 60.50d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with CH (1-5) [jpcb2004]
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.92d0
      nbclocal(2) = 40.00d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with C (1-6) [jpcb2004]
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 4.56d0
      nbclocal(2) = 10.00d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with Na (1-7) [jacs2004]
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.4d0
      nbclocal(2) = 23.0d0
      call twh_set_nbcoeff(ij,nbclocal) 

c     Methane (CH4) united-atom [jpcb2004]
      type = 2
      masslocal = 12.001d0 + 4.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'ion'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of CH4 with CH4 (2-2) [jpcb2004]
      j = 2
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.72d0
      nbclocal(2) = 158.50d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH4 with CH3 (2-3) [jpcb2004]
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.74d0
      nbclocal(2) = 130.84d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH4 with CH2 (2-4) [jpcb2004]
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.84d0
      nbclocal(2) = 94.21d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH4 with CH (2-5) [jpcb2004]
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 4.17d0
      nbclocal(2) = 51.91d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH4 with C (2-6) [jpcb2004]
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 4.87d0
      nbclocal(2) = 11.26d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH4 with Na (2-7) [jacs2004]
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.72d0
      nbclocal(2) = 582.17d0
      call twh_set_nbcoeff(ij,nbclocal) 

c     Methyl (CH3) united-atom [jpcb2004]
      type = 3
      masslocal = 12.001d0 + 3.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3'
      cbnlocal = 'Csp3'
      canlocal = 'Csp3'
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of CH3 with CH3 (3-3) [jpcb2004]
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.76d0
      nbclocal(2) = 108.00d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH3 with CH2 (3-4) [jpcb2004]
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.86d0
      nbclocal(2) = 77.77d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH3 with CH (3-5) [jpcb2004]
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 4.19d0
      nbclocal(2) = 42.85d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH3 with C (3-6) [jpcb2004]
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 4.90d0
      nbclocal(2) = 9.30d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH3 with C (3-7) [jacs2004]
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.65d0
      nbclocal(2) = 443.73d0
      call twh_set_nbcoeff(ij,nbclocal) 

c     methylene (CH2) united-atom [jpcb2004]
      type = 4
      masslocal = 12.001d0 + 2.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2'
      cbnlocal = 'Csp3'
      canlocal = 'Csp3'
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of CH2 with CH2 (4-4) [jpcb2004]
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.96d0
      nbclocal(2) = 56.00d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH2 with CH (4-5) [jpcb2004]
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 4.30d0
      nbclocal(2) = 30.85d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH2 with C (4-6) [jpcb2004]
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 5.03d0
      nbclocal(2) = 6.69d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH2 with Na (4-7) [jacs2004]
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.95d0
      nbclocal(2) = 310.0d0
      call twh_set_nbcoeff(ij,nbclocal) 

c     methine (CH) united-atom [jpcb2004]
      type = 5
      masslocal = 12.001d0 + 1.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH'
      cbnlocal = 'Csp3'
      canlocal = 'Csp3'
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of CH with CH (5-5) [jpcb2004]
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 4.67d0
      nbclocal(2) = 17.00d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of CH with C (5-6) [jpcb2004]
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 5.46d0
      nbclocal(2) = 3.69d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- CH with Na (5-7) not yet available

c     quaternary carbon (C) [jpcb2004]
      type = 6
      masslocal = 12.001d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C'
      cbnlocal = 'Csp3'
      canlocal = 'Csp3'
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of C with C (6-6) [jpcb2004]
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 6.38d0
      nbclocal(2) = 0.80d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- C with Na (6-7) not yet available

c     sodium ion [jacs2004]
      type = 7
      masslocal = 22.990d0
      elelocal =  'Na'
      bplocal = 'ion'
      call twh_qbasevalue(GLB_SET,type,1.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Na'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of Na with Na (7-7) [jacs2004]
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.16d0
      nbclocal(2) = 124.4d0
      call twh_set_nbcoeff(ij,nbclocal) 

c     aluminum in zeolite framework [jacs2004]
      type = 8
      masslocal = 26.982d0
      elelocal =  'Al'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,-1.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Al'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- does not interact via lennard-jones, only coulombic
c     --- no interacts set with any other atom [jacs2004]

c     silicon in zeolite framework [jacs2004]
      type = 9
      masslocal = 28.086d0
      elelocal =  'Si'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Si'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- does not interact via lennard-jones, only coulombic
c     --- no interacts set with any other atom [jacs2004]

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond lengths - simple harmonic                                 *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.54d0
      vclocal(1) = (0.5d0)*(96500d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- any bond between sp3 carbons
      isame = isame + 1
      vnlocal(isame,1) = 'Csp3'
      vnlocal(isame,2) = 'Csp3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     ******************************************************************
c     * more bond lengths - no interaction                             *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_NONE)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- Al-O bonds in a zeolite
      isame = isame + 1
      vnlocal(isame,1) = 'Al'
      vnlocal(isame,2) = 'O'
c     --- Si-O bonds in a zeolite
      isame = isame + 1
      vnlocal(isame,1) = 'Si'
      vnlocal(isame,2) = 'O'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * angles - harmonic Cosine                                       *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 114.0d0
      bclocal(1) = (0.5d0)*(62500.0d0)
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     all alkanes [jpcb2004]
      isame = isame + 1
      bnlocal(isame,1) = 'Csp3'
      bnlocal(isame,2) = 'Csp3'
      bnlocal(isame,3) = 'Csp3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     ******************************************************************
c     * more angles - no interaction                                   *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_NONE)
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     Si - O - Si in zeolite
      isame = isame + 1
      bnlocal(isame,1) = 'Si'
      bnlocal(isame,2) = 'O'
      bnlocal(isame,3) = 'Si'
c     Si - O - Al in zeolite
      isame = isame + 1
      bnlocal(isame,1) = 'Si'
      bnlocal(isame,2) = 'O'
      bnlocal(isame,3) = 'Al'
c     O - Si - O  in zeolite
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'Si'
      bnlocal(isame,3) = 'O'
c     O - Al - O  in zeolite
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'Al'
      bnlocal(isame,3) = 'O'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * torsions - cosine power series                                 *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,10)
      call twh_ntorloop(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 1204.654d0
      tclocal(1) = 1947.740d0
      tclocal(2) = -357.845d0
      tclocal(3) = -1944.666d0
      tclocal(4) = 715.690d0
      tclocal(5) = -1565.572d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     type 1, example n-butane
c     C_x - CH2 - CH2 - C_x [jpcb2004]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2'
      tnlocal(isame,3) = 'CH2'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,10)
      call twh_ntorloop(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 1293.324d0
      tclocal(1) = 3879.849d0
      tclocal(2) = 0.0d0
      tclocal(3) = -5173.163d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     type 3, example 2,2-dimethylbutane
c     C_x - C - CH2 - C_x [jpcb2004]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CH2'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,10)
      call twh_ntorloop(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 2045.657d0
      tclocal(1) = 6136.797d0
      tclocal(2) = 0.0d0
      tclocal(3) = -8182.447d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     type 4, example 2,2,3,3-tetramethylbutane
c     C_x - C - C - C_x [jpcb2004]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     ******************************************************************
c     * more torsions, these are stolen from the TraPPE force field    *
c     * and are used in place of the rather bizarre Dubb2004 torsions  *
c     * that involve a hydrogen represented on the united-atom of CH   *
c     *                                                                *
c     * replacing with TraPPE torsions as suggested by T.J.H. Vlugt    *
c     * [personal communication M.G. Martin with T.J.H. Vlught 12/2004 *
c     ******************************************************************

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
c     X - CH - CH2 - X taken from TraPPE
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'CH2'
      tnlocal(isame,4) = 'Cwild'
c     X - CH - CH - X taken from TraPPE
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH'
      tnlocal(isame,3) = 'CH'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     ******************************************************************
c     * more torsions - null interactions                              *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,8)
      call twh_ntorloop(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     X - Al - O - X in zeolite
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'Al'
      tnlocal(isame,3) = 'O'
      tnlocal(isame,4) = 'wild'
c     X - Si - O - X bonds in zeolite
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'Si'
      tnlocal(isame,3) = 'O'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond increments                                                *
c     ******************************************************************

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.5125d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
      qbinlocal(1) = 'Si'
      qbinlocal(2) = 'O'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.6875d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
      qbinlocal(1) = 'Al'
      qbinlocal(2) = 'O'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Csp3'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
