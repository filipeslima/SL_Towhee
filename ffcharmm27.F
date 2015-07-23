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
      subroutine twh_ffcharmm27(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the charmm27     *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-20-2002 by M.G. Martin  *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
c     --- local variables
      character*(FFNAMELEN) ffname

      ffname = 'Charmm27'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      call twh_ffcharmm27_nonbond(lfailure,ffname,ntypes)
      if ( lfailure ) return
      call twh_ffcharmm27_bond(lfailure,ffname,nbondtypes)
      if ( lfailure ) return
      call twh_ffcharmm27_angle(lfailure,ffname,nangletypes)
      if ( lfailure ) return
      call twh_ffcharmm27_dihedral(lfailure,ffname,ndihedtypes)
      if ( lfailure ) return
      call twh_ffcharmm27_improper(lfailure,ffname,nimprotypes)
      if ( lfailure ) return

      return
      end


      subroutine twh_ffcharmm27_nonbond(lfailure,ffname,ntypes)
c     ******************************************************************
c     * sets up the nonbonded portion of the charmm27 forcefield       *
c     *                                                                *
c     * split from the original ffcharmm27 03-05-2006 by M.G. Martin   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer ntypes
c     --- local variables
      integer type,ii
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

c     **************************************************************
c     * Charmm27 force field                                       *
c     **************************************************************
c     NOTE: the charmm27 parameters were obtained from Alex 
c     MacKerell's web page
c   http://www.pharmacy.umaryland.edu/faculty/amackere/force_fields.htm
c     the appropriate references can be found at the top of those
c     files.

c     --- initialize type count
      type = 0
c     *** Calcium ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.367d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.120d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 40.08000d0 
      bplocal = 'ion'
      elelocal = 'Ca'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CAL) calcium ion [par_all27_na.inp]
      nbnlocal ='CAL'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Carbon ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.11d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (C) carbonyl C, peptide backbone [par_all27_prot_na.inp]
      nbnlocal ='C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9924d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.07d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CA) aromatic C [par_all27_prot_na.inp]
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
c     (CC)carbonyl C,asn,asp,gln,glu,cter,ct2 [par_all27_prot_na.inp]
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
c     (CD) carbonyl C, pres aspp,glup,ct1 [par_all27_prot_na.inp]
      nbnlocal ='CD'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.09d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.068d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CE1) for propene [par_all27_prot_lipid.inp]
      nbnlocal ='CEL1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.08d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.064d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      bplocal = 'sp2'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CE2) for alkene; H2C=CR [par_all27_prot_lipid.inp]
      nbnlocal ='CE2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.09d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.068d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CEL1) alkene [par_all27_lipid.inp]
      nbnlocal ='CEL1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.08d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.064d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      bplocal = 'sp2'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CEL2) for alkene; H2C=CR [par_all27_lipid.inp]
      nbnlocal ='CEL2'
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
c     (CL) methyl acetate update [par_all27_lipid.inp]
      nbnlocal ='CL'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.100d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.11d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CM) heme CO carbon [par_all27_prot_lipid.inp]
      nbnlocal ='CM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.10d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN1) [par_all27_na.inp]
      nbnlocal ='CN1'
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
c     (CN1A) [par_all27_na.inp]
      nbnlocal ='CN1A'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.10d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN1T) [par_all27_na.inp]
      nbnlocal ='CN1T'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.10d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN2) [par_all27_na.inp]
      nbnlocal ='CN2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.09d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN3) [par_all27_na.inp]
      nbnlocal ='CN3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.800d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.18d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN3A) [par_all27_na.inp]
      nbnlocal ='CN3A'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.800d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.18d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN3B) [par_all27_na.inp]
      nbnlocal ='CN3B'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.800d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.18d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN3C) [par_all27_na.inp]
      nbnlocal ='CN3C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.09d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN3D) [par_all27_na.inp]
      nbnlocal ='CN3D'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.09d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN3T) [par_all27_na.inp]
      nbnlocal ='CN3T'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.075d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN4) [par_all27_na.inp]
      nbnlocal ='CN4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.075d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN5) [par_all27_na.inp]
      nbnlocal ='CN5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.075d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN5G) [par_all27_na.inp]
      nbnlocal ='CN5G'
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
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN7) [par_all27_na.inp]
      nbnlocal ='CN7'
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
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN7B) [par_all27_na.inp]
      nbnlocal ='CN7B'
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
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN7C) [par_all27_na.inp]
      nbnlocal ='CN7C'
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
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN7D) [par_all27_na.inp]
      nbnlocal ='CN7D'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.010d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0560d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.01d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
c     (CN8) [par_all27_na.inp]
      nbnlocal ='CN8'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.010d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0560d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.01d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
c     (CN8B) [par_all27_na.inp]
      nbnlocal ='CN8B'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.040d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0780d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.01d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
c     (CN9) [par_all27_na.inp]
      nbnlocal ='CN9'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9924d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.07d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CNA) [par_all27_na.inp]
      nbnlocal ='CNA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CNA2) [par_all27_na.inp]
      nbnlocal ='CNA2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.09d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.068d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CNE1) [par_all27_na.inp]
      nbnlocal ='CNE1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.08d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.064d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      bplocal = 'sp2'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CNE2) [par_all27_na.inp]
      nbnlocal ='CNE2'
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
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CP1) tetrahedral C (proline CA) [par_all27_prot_lipid.inp]
      nbnlocal ='CP1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.175d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.055d0)*(CON_KCALTOK)
      nbclocal(3) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.010d0)*(CON_KCALTOK)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CP2) tetrahedral C (proline CB/CG) [par_all27_prot_lipid.inp]
      nbnlocal ='CP2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.175d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.055d0)*(CON_KCALTOK)
      nbclocal(3) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.010d0)*(CON_KCALTOK)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CP3) tetrahedral C (proline CB/CG) [par_all27_prot_lipid.inp]
      nbnlocal ='CP3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.800d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.09d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPA) heme alpha-C [par_all27_prot_lipid.inp]
      nbnlocal ='CPA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.800d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.09d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPB) heme beta-C [par_all27_prot_lipid.inp]
      nbnlocal ='CPB'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.800d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPH1) his CG and CD2 carbons [par_all27_prot_lipid.inp]
      nbnlocal ='CPH1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.800d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPH2) his CE1 carbon [par_all27_prot_lipid.inp]
      nbnlocal ='CPH2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.800d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.09d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPM) heme meso-C [par_all27_prot_lipid.inp]
      nbnlocal ='CPM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.800d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.090d0)*(CON_KCALTOK)
      nbclocal(3) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.090d0)*(CON_KCALTOK)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPT) trp C between rings [par_all27_prot_lipid.inp]
      nbnlocal ='CPT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.200d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.11d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CS) thiolate carbon [par_all27_prot_lipid.inp]
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
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT1) alkane [par_all27_prot_lipid.inp]
      nbnlocal ='CT1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.175d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.055d0)*(CON_KCALTOK)
      nbclocal(3) = (1.900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.010d0)*(CON_KCALTOK)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT2) aliphatic sp3 C for CH2 [par_all27_prot_lipid.inp]
      nbnlocal ='CT2'
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
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT3) aliphatic sp3 C for CH3 [par_all27_prot_lipid.inp]
      nbnlocal ='CT3'
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
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CTL1) alkane [par_all27_lipid.inp]
      nbnlocal ='CTL1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.010d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0560d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.01d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
c     (CTL2) alkane [par_all27_lipid.inp]
      nbnlocal ='CTL2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.040d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0780d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.01d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
c     (CTL3) alkane [par_all27_lipid.inp]
      nbnlocal ='CTL3'
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
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CTL5) tetramethylammonium [par_all27_lipid.inp]
      nbnlocal ='CTL5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.9924d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.07d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.000d0
      bplocal = 'sp3'
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CY) TRP C in pyrrole ring [par_all27_prot_lipid.inp]
      nbnlocal ='CY'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Cesium ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.100d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.190d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 132.90d0 
      bplocal = 'ion'
      elelocal = 'Cs'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CES) cesium ion [par_all27_na.inp]
      nbnlocal ='CES'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Chlorine *** 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.270d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.150d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 35.45000d0 
      bplocal = 'ion'
      elelocal = 'Cl'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      call twh_qbasevalue(GLB_SET,type,-1.0d0)
c     (CLA) chloride ion [par_all27_na.inp]
      nbnlocal ='CLA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Fluorine ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.700d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.120d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 18.9984d0
      bplocal = 'sp3'
      elelocal = ' F'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (FA1) aromatic fluorine [par_all27_na.inp]
      nbnlocal ='FA1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.700d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.090d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 18.9984d0
      bplocal = 'sp3'
      elelocal = ' F'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (FN1) Fluorine for sugar derivatives [par_all27_na.inp]
      nbnlocal ='FN1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Hydrogen ***

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
c     (H) polar H [par_all27_prot_lipid.inp]
      nbnlocal ='H'
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
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HA) nonpolar H [par_all27_prot_lipid.inp]
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
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HAL1) [par_all27_lipid.inp]
      nbnlocal = 'HAL1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3400d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.028d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
c     (HAL2) [par_all27_lipid.inp]
      nbnlocal = 'HAL2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3400d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.024d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
c     (HAL3) [par_all27_lipid.inp]
      nbnlocal = 'HAL3'
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
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HB) backbone H [par_all27_prot_lipid.inp]
      nbnlocal ='HB'
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
c     (HC) N-ter H [par_all27_prot_lipid.inp]
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
c     (HCL) ethnaolamine [par_all27_lipid.inp]
      nbnlocal ='HCL'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.25d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.031d0)*(CON_KCALTOK)
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
c     (HE1) alkene; RHC=CR [par_all27_prot_lipid.inp]
      nbnlocal ='HE1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.26d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0261d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HE2) alkene; H2C=CR [par_all27_prot_lipid.inp]
      nbnlocal ='HE2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.25d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.031d0)*(CON_KCALTOK)
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
c     (HEL1) alkene [par_all27_lipid.inp]
      nbnlocal ='HEL1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.26d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0261d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HEL2) alkene [par_all27_lipid.inp]
      nbnlocal ='HEL2'
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
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HL) polar H on NC4+ [par_all27_lipid.inp]
      nbnlocal ='HL'
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
c     (HN1) [par_all27_na.inp]
      nbnlocal ='HN1'
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
c     (HN2) [par_all27_na.inp]
      nbnlocal ='HN2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.1000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
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
c     (HN3) [par_all27_na.inp]
      nbnlocal ='HN3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.9000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0460d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HN3B) [par_all27_na.inp]
      nbnlocal ='HN3B'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3582d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HN3C) [par_all27_na.inp]
      nbnlocal ='HN3C'
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
c     (HN4) [par_all27_na.inp]
      nbnlocal ='HN4'
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
c     (HN5) [par_all27_na.inp]
      nbnlocal ='HN5'
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
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HN6) [par_all27_na.inp]
      nbnlocal ='HN6'
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
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HN7) [par_all27_na.inp]
      nbnlocal ='HN7'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3400d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.028d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
c     (HN8) [par_all27_na.inp]
      nbnlocal ='HN8'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3400d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.024d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
c     (HN9) [par_all27_na.inp]
      nbnlocal ='HN9'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.25d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.031d0)*(CON_KCALTOK)
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
c     (HNE1) [par_all27_na.inp]
      nbnlocal ='HNE1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.26d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0260d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HNE2) [par_all27_na.inp]
      nbnlocal ='HNE2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3582d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HNP) [par_all27_na.inp]
      nbnlocal ='HNP'
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
c     (HOL) [par_all27_lipid.inp]
      nbnlocal ='HOL'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.3582d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HP) aromatic H [par_all27_prot_lipid.inp]
      nbnlocal ='HP'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.9000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0460d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HR1) his he1, (+) his HG,HD2 [par_all27_prot_lipid.inp]
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
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HR2) (+) his HE1 [par_all27_prot_lipid.inp]
      nbnlocal ='HR2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.468d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.00780d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HR3) neutral his HG, HD2 [par_all27_prot_lipid.inp]
      nbnlocal ='HR3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.450d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      bplocal = 's'
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HS) thiol hydrogen [par_all27_prot_lipid.inp]
      nbnlocal ='HS'
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
c     (HT) water hydrogen [par_all27_lipid.inp]
      nbnlocal ='HT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Iron ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.6500d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 132.90d0 
      bplocal = 'ion'
      elelocal = 'Fe'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (FE) iron in heme [par_all27_prot_lipid.inp]
      nbnlocal ='FE'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Magnesium *** 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.185d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.015d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 24.30500d0 
      bplocal = 'ion'
      elelocal = 'Mg'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      call twh_qbasevalue(GLB_SET,type,2.0d0)
c     (MG) magnesium ion [par_all27_na.inp]
      nbnlocal ='MG'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Nitrogen ***
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.85d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = (1.85d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.0001d0)*(CON_KCALTOK)
      masslocal = 14.007d0
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (N) proline N [par_all27_prot_lipid.inp]
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NC2) guanidinium nitroogen [par_all27_prot_lipid.inp]
      nbnlocal ='NC2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.85d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2000d0)*(CON_KCALTOK)
      nbclocal(3) = (1.55d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.2000d0)*(CON_KCALTOK)
      masslocal = 14.007d0
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NH1) peptide nitrogen [par_all27_prot_lipid.inp]
      nbnlocal ='NH1'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NH2) amide nitrogen [par_all27_prot_lipid.inp]
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NH3) ammonium nitrogen [par_all27_prot_lipid.inp]
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NH3L) ethanolamine [par_all27_lipid.inp]
      nbnlocal ='NH3L'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN1) [par_all27_na.inp]
      nbnlocal ='NN1'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN1C) [par_all27_na.inp]
      nbnlocal ='NN1C'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN2) [par_all27_na.inp]
      nbnlocal ='NN2'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN2B) [par_all27_na.inp]
      nbnlocal ='NN2B'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN2C) [par_all27_na.inp]
      nbnlocal ='NN2C'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN2G) [par_all27_na.inp]
      nbnlocal ='NN2G'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN2U) [par_all27_na.inp]
      nbnlocal ='NN2U'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN3) [par_all27_na.inp]
      nbnlocal ='NN3'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN3A) [par_all27_na.inp]
      nbnlocal ='NN3A'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN3G) [par_all27_na.inp]
      nbnlocal ='NN3G'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN3I) [par_all27_na.inp]
      nbnlocal ='NN3I'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN4) [par_all27_na.inp]
      nbnlocal ='NN4'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN5) [par_all27_na.inp]
      nbnlocal ='NN5'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NN6) [par_all27_na.inp]
      nbnlocal ='NN6'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NP) Proline ring NH2+ (N-terminal) [par_all27_prot_lipid.inp]
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NPH) heme pyrrole N [par_all27_prot_lipid.inp]
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NR1) neutral his protonated ring nitrogen 
c      --- [par_all27_prot_lipid.inp]
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NR2) neutral his unprotonated ring nitrogen 
c      --- [par_all27_prot_lipid.inp]
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NR3) charged his ring nitrogen [par_all27_prot_lipid.inp]
      nbnlocal ='NR3'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NTL) as all other nitrogens [par_all27_lipid.inp]
      nbnlocal ='NTL'
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
      bplocal = 'sp3'
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NY) TRP N in pyrrole ring [par_all27_prot_lipid.inp]
      nbnlocal ='NY'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Oxygen ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1200d0)*(CON_KCALTOK)
      nbclocal(3) = (1.400d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.1200d0)*(CON_KCALTOK)
      masslocal = 15.999d0
      bplocal = 's'
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (O) carbonyl oxygen [par_all27_prot_lipid.inp]
      nbnlocal ='O'
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
c     (O2L) [par_all27_lipid.inp]
      nbnlocal ='O2L'
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
      bplocal = 's'
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OB) carbonyl oxygen in acetic acid [par_all27_prot_lipid.inp]
      nbnlocal ='OB'
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
      bplocal = 's'
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OBL) [par_all27_lipid.inp]
      nbnlocal ='OBL'
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
c     (OC) carboxylate oxygen [par_all27_prot_lipid.inp]
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
c     (OCL) [par_all27_lipid.inp]
      nbnlocal ='OCL'
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
c     (OH1) hydroxyl oxygen [par_all27_prot_lipid.inp]
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
c     (OHL) [par_all27_lipid.inp]
      nbnlocal ='OHL'
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
c     (OM) heme CO/O2 oxygen [par_all27_prot_lipid.inp]
      nbnlocal ='OM'
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
c     (ON1) [par_all27_na.inp]
      nbnlocal ='ON1'
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
c     (ON1C) [par_all27_na.inp]
      nbnlocal ='ON1C'
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
c     (ON2) [par_all27_na.inp]
      nbnlocal ='ON2'
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
c     (ON3) [par_all27_na.inp]
      nbnlocal ='ON3'
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
c     (ON4) [par_all27_na.inp]
      nbnlocal ='ON4'
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
c     (ON5) [par_all27_na.inp]
      nbnlocal ='ON5'
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
c     (ON6) [par_all27_na.inp]
      nbnlocal ='ON6'
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
c     (ON6B) [par_all27_na.inp]
      nbnlocal ='ON6B'
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
c     (OS) ester oxygen [par_all27_prot_lipid.inp]
      nbnlocal ='OS'
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
c     (OSL) [par_all27_lipid.inp]
      nbnlocal ='OSL'
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
c     (OT) TIP3P water oxygen [par_all27_lipid.inp]
      nbnlocal ='OT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Phosphorus ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.15d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.585d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 30.974d0
      elelocal = ' P'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
c     (P)  [par_all27_na.inp]
      nbnlocal ='P'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.15d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.585d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 30.974d0
      elelocal = ' P'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
c     (P2)  [par_all27_na.inp]
      nbnlocal ='P2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.15d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.585d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 30.974d0
      elelocal = ' P'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
c     (PL)  [par_all27_lipid.inp]
      nbnlocal ='PL'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Potasium ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.76375d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.087d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 39.10200d0 
      bplocal = 'ion'
      elelocal = ' K'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      call twh_qbasevalue(GLB_SET,type,1.0d0)
c     (POT) potassium ion [par_all27_na.inp]
      nbnlocal ='POT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Sodium ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.36375d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0469d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 22.98977d0 
      bplocal = 'ion'
      elelocal = 'Na'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      call twh_qbasevalue(GLB_SET,type,1.0d0)
c     (SOD) sodium ion [par_all27_na.inp]
      nbnlocal ='SOD'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Sulfur ***

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.00d0)*(2.0d0)/(2.0d0**(1.0/6.0))
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
c     (S) sulphur [par_all27_prot_lipid.inp]
      nbnlocal ='S'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.1d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.47d0)*(CON_KCALTOK)
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
c     (SL) methylsulfate [par_all27_lipid.inp]
      nbnlocal ='SL'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.975d0)*(2.0d0)/(2.0d0**(1.0/6.0))
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
c     (SM) sulfur C-S-S-C type [par_all27_prot_lipid.inp]
      nbnlocal ='SM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.200d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.470d0)*(CON_KCALTOK)
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
c     (SS) thiolate sulfur [par_all27_prot_lipid.inp]
      nbnlocal ='SS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffcharmm27_bond(lfailure,ffname,nbondtypes)
c     ******************************************************************
c     * sets up the bond potentials for Charmm27                       *
c     *                                                                *
c     * split from original ffcharmm27 03-05-2006 by M.G. Martin       *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nbondtypes
c     --- local variables
      integer type,isame
      integer nvslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)

c     --- initialize type count
      type = 0

c     ******************************************************************
c     * Charmm27 bond potentials                                       *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (200.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL3 - CL methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL3'
      vnlocal(isame,2) = 'CL'
c     --- CTL2 - CL methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL2'
      vnlocal(isame,2) = 'CL'
c     --- CTL1 - CL methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL1'
      vnlocal(isame,2) = 'CL'
c     --- CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL1'
      vnlocal(isame,2) = 'CC'
c     --- CT1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL1'
      vnlocal(isame,2) = 'CD'
c     --- CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'CC'
c     --- CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'CD'
c     --- CT3 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'CC'
c     --- CT3 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'CD'
c     --- CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT1'
      vnlocal(isame,2) = 'CC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.500d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL1 - CTL1 alkanes [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL1'
      vnlocal(isame,2) = 'CTL1'
c     --- CTL - CTL alkanes [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT1'
      vnlocal(isame,2) = 'CT1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5380d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL1 - CTL2 alkanes [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL1'
      vnlocal(isame,2) = 'CTL2'
c     --- CTL1 - CTL3 alkanes [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL1'
      vnlocal(isame,2) = 'CTL3'
c     --- CT3 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'CT1'
c     --- CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'CT1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5280d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL2 - CTL3 alkanes [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL2'
      vnlocal(isame,2) = 'CTL3'
c     --- CT3 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'CT2'
c     --- CN8 - CN9 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8'
      vnlocal(isame,2) = 'CN9'
c     --- CN8 - CN8 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8'
      vnlocal(isame,2) = 'CN8'
c     ---- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5300d0
      vclocal(1) = (222.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL2 - CTL2 alkanes [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL2'
      vnlocal(isame,2) = 'CTL2'
c     --- CTL3 - CTL3 alkanes [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL3'
      vnlocal(isame,2) = 'CTL3'
c     --- CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'CT2'
c     --- CT3 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'CT3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1110d0
      vclocal(1) = (309.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL1 - HAL1 alkanes [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL1'
      vnlocal(isame,2) = 'HAL1'
c     --- CTL2 - HAL2 alkanes [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL2'
      vnlocal(isame,2) = 'HAL2'
c     --- HA - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CP2'
c     --- HA - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CP3'
c     --- HA - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CT2'
c     --- CN7 - HN7 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'HN7'
c     --- CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'HN7'
c     --- CN7C - HN7 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7C'
      vnlocal(isame,2) = 'HN7'
c     --- CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8'
      vnlocal(isame,2) = 'HN8'
c     --- CN8B - HN8 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8B'
      vnlocal(isame,2) = 'HN8'
c     --- HA - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CT1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1110d0
      vclocal(1) = (322.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL3 - HAL3 alkanes [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL3'
      vnlocal(isame,2) = 'HAL3'
c     --- HA3 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CT3'
c     --- CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN9'
      vnlocal(isame,2) = 'HN9'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.220d0
      vclocal(1) = (750.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OBL - CL methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OBL'
      vnlocal(isame,2) = 'CL'
c     --- OB - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OB'
      vnlocal(isame,2) = 'CC'
c     --- OB - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OB'
      vnlocal(isame,2) = 'CD'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.2600d0
      vclocal(1) = (525.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OCL - CL acetate, protein [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OCL'
      vnlocal(isame,2) = 'CL'
c     --- OC - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OC'
      vnlocal(isame,2) = 'CA'
c     --- OC - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OC'
      vnlocal(isame,2) = 'CC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4000d0
      vclocal(1) = (230.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OHL - CL methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OHL'
      vnlocal(isame,2) = 'CL'
c     --- OH1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OH1'
      vnlocal(isame,2) = 'CD'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4200d0
      vclocal(1) = (428.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OHL - CTL1 glycerol [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OHL'
      vnlocal(isame,2) = 'CTL1'
c     --- OHL - CTL2 glycerol [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OHL'
      vnlocal(isame,2) = 'CTL2'
c     --- OHL - CTL3 glycerol [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OHL'
      vnlocal(isame,2) = 'CTL3'
c     --- OH1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OH1'
      vnlocal(isame,2) = 'CT1'
c     --- OH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OH1'
      vnlocal(isame,2) = 'CT2'
c     --- OH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OH1'
      vnlocal(isame,2) = 'CT3'
c     --- CN7 - ON5 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'ON5'
c     --- CN8 - ON5 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8'
      vnlocal(isame,2) = 'ON5'
c     --- CN8B - ON5 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8B'
      vnlocal(isame,2) = 'ON5'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9600d0
      vclocal(1) = (545.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HOL - OHL acetic acid [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HOL'
      vnlocal(isame,2) = 'OHL'
c     --- OH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OH1'
      vnlocal(isame,2) = 'H'
c     --- ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'ON5'
      vnlocal(isame,2) = 'HN5'
c     --- HN4 - ON4 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HN4'
      vnlocal(isame,2) = 'ON4'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3340d0
      vclocal(1) = (150.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OSL - CL methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OSL'
      vnlocal(isame,2) = 'CL'
c     --- OS - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OS'
      vnlocal(isame,2) = 'CD'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4300d0
      vclocal(1) = (340.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL3 - OSL phosphate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL3'
      vnlocal(isame,2) = 'OSL'
c     --- CTL2 - OSL phosphate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL2'
      vnlocal(isame,2) = 'OSL'
c     --- CTL1 - OSL phosphate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL1'
      vnlocal(isame,2) = 'OSL'
c     --- OS - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OS'
      vnlocal(isame,2) = 'CT3'
c     --- CN9 - ON2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN9'
      vnlocal(isame,2) = 'ON2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9572d0
      vclocal(1) = (450.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OT - HT from TIPS3P geometry [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OT'
      vnlocal(isame,2) = 'HT'
c     --- set nvsame
      nvslocal = 1
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.60d0
      vclocal(1) = (270.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OSL - PL phosphate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OSL'
      vnlocal(isame,2) = 'PL'
c     --- ON2 - P [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'ON2'
      vnlocal(isame,2) = 'P'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.48d0
      vclocal(1) = (580.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- O2L - PL phosphate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'O2L'
      vnlocal(isame,2) = 'PL'
c     --- ON3 - P [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'ON3'
      vnlocal(isame,2) = 'P'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.59d0
      vclocal(1) = (237.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OHL - PL phosphate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OHL'
      vnlocal(isame,2) = 'PL'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.04d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH3L - HCL ethanolamine [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH3L'
      vnlocal(isame,2) = 'HCL'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.51d0
      vclocal(1) = (261.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH3L - CTL2 ethanolamine [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH3L'
      vnlocal(isame,2) = 'CTL2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.51d0
      vclocal(1) = (215.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NTL - CTL2 tetramethylammonium [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NTL'
      vnlocal(isame,2) = 'CTL2'
c     --- NTL - CTL5 tetramethylammonium [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NTL'
      vnlocal(isame,2) = 'CTL5'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (300.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL5 - HL tetramethylammonium [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL5'
      vnlocal(isame,2) = 'HL'
c     --- CTL2 - HL tetramethylammonium [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CTL2'
      vnlocal(isame,2) = 'HL'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.448d0
      vclocal(1) = (540.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- SL - O2L methylsulfate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'SL'
      vnlocal(isame,2) = 'O2L'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.575d0
      vclocal(1) = (250.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- SL - OSL methylsulfate [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'SL'
      vnlocal(isame,2) = 'OSL'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.330d0
      vclocal(1) = (510.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL2 - CEL2 ethene [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CEL2'
      vnlocal(isame,2) = 'CEL2'
c     --- CE2 - CE2 for ethene [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CE2'
      vnlocal(isame,2) = 'CE2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.100d0
      vclocal(1) = (365.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL2 - CEL2 propene [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HEL2'
      vnlocal(isame,2) = 'CEL2'
c     --- HE2 - CE2 for ethene [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HE2'
      vnlocal(isame,2) = 'CE2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.504d0
      vclocal(1) = (383.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL1 - CTL3 butene [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CEL1'
      vnlocal(isame,2) = 'CTL3'
c     --- CE1 - CT3 for butene [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CE1'
      vnlocal(isame,2) = 'CT3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.342d0
      vclocal(1) = (500.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL1 - CEL2 propene [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CEL1'
      vnlocal(isame,2) = 'CEL2'
c     --- CE1 - CE2 for propene [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CE1'
      vnlocal(isame,2) = 'CE2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.100d0
      vclocal(1) = (360.5d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL1 - CEL1 propene [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HEL1'
      vnlocal(isame,2) = 'CEL1'
c     --- HE1 - CE1 for propene [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HE1'
      vnlocal(isame,2) = 'CE1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.502d0
      vclocal(1) = (365.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL1 - CTL2 butene [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CEL1'
      vnlocal(isame,2) = 'CTL2'
c     --- CE1 - CT2 for butene; from propene [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CE1'
      vnlocal(isame,2) = 'CT2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.340d0
      vclocal(1) = (440.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL1 - CEL1 butene [par_all27_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CEL1'
      vnlocal(isame,2) = 'CEL1'
c     --- CE1 - CE1 for butene; from propene [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CE1'
      vnlocal(isame,2) = 'CE1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     *** Protein section (for protein bond groups not listed above) 
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (600.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- C - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.375d0
      vclocal(1) = (305.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CA'
      vnlocal(isame,2) = 'CA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.490d0
      vclocal(1) = (250.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CP1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CP1'
      vnlocal(isame,2) = 'C'
c     --- CP1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CP1'
      vnlocal(isame,2) = 'CC'
c     --- CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT1'
      vnlocal(isame,2) = 'C'
c     --- CT2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'C'
c     --- CT3 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.490d0
      vclocal(1) = (200.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CP1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CP1'
      vnlocal(isame,2) = 'CD'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.527d0
      vclocal(1) = (222.5d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CP2 - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CP2'
      vnlocal(isame,2) = 'CP1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.537d0
      vclocal(1) = (222.5d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CP2 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CP2'
      vnlocal(isame,2) = 'CP2'
c     --- CP3 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CP3'
      vnlocal(isame,2) = 'CP2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.380d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPB - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CPB'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame 
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4432d0
      vclocal(1) = (299.8d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPB - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CPB'
      vnlocal(isame,2) = 'CPA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3464d0
      vclocal(1) = (340.7d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPB - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CPB'
      vnlocal(isame,2) = 'CPB'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3600d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CPH1'
      vnlocal(isame,2) = 'CPH1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3716d0
      vclocal(1) = (360.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPM - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CPM'
      vnlocal(isame,2) = 'CPA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3680d0
      vclocal(1) = (305.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CPT'
      vnlocal(isame,2) = 'CA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4000d0
      vclocal(1) = (360.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPT - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CPT'
      vnlocal(isame,2) = 'CPT'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4900d0
      vclocal(1) = (230.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'CA'
c     --- CT2 - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'CPB'
c     --- CT3 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'CA'
c     --- CT3 - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'CPB'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5000d0
      vclocal(1) = (229.63d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'CPH1'
c     --- CT3 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'CPH1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5310d0
      vclocal(1) = (190.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CT3 - CS [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'CS'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3650d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CY'
      vnlocal(isame,2) = 'CA'
c     --- CN4 - NN2G [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN4'
      vnlocal(isame,2) = 'NN2G'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4400d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CY'
      vnlocal(isame,2) = 'CPT'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5100d0
      vclocal(1) = (230.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CY - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CY'
      vnlocal(isame,2) = 'CT2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1100d0
      vclocal(1) = (330.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- H - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'H'
      vnlocal(isame,2) = 'CD'
c     set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1000d0
      vclocal(1) = (330.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0830d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1000d0
      vclocal(1) = (317.13d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0900d0
      vclocal(1) = (367.6d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - CPM [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CPM'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1110d0
      vclocal(1) = (300.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - CS [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CS'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0800d0
      vclocal(1) = (330.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HA'
      vnlocal(isame,2) = 'CY'
c     --- HB - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HB'
      vnlocal(isame,2) = 'CP1'
c     --- HB - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HB'
      vnlocal(isame,2) = 'CT1'
c     --- HB - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HB'
      vnlocal(isame,2) = 'CT2'
c     --- HB - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HB'
      vnlocal(isame,2) = 'CT3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0800d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HP - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HP'
      vnlocal(isame,2) = 'CA'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0800d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HP - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HP'
      vnlocal(isame,2) = 'CY'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0830d0
      vclocal(1) = (375.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HR1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HR1'
      vnlocal(isame,2) = 'CPH1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0900d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HR1 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HR1'
      vnlocal(isame,2) = 'CPH2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0700d0
      vclocal(1) = (333.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HR2 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HR2'
      vnlocal(isame,2) = 'CPH2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0830d0
      vclocal(1) = (365.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HR3 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HR3'
      vnlocal(isame,2) = 'CPH1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3000d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- N - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'N'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4340d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- N - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'N'
      vnlocal(isame,2) = 'CP1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4550d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- N - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'N'
      vnlocal(isame,2) = 'CP3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3650d0
      vclocal(1) = (463.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NC2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NC2'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4900d0
      vclocal(1) = (261.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NC2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NC2'
      vnlocal(isame,2) = 'CT2'
c     --- NC2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NC2'
      vnlocal(isame,2) = 'CT3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (455.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NC2 - HC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NC2'
      vnlocal(isame,2) = 'HC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3450d0
      vclocal(1) = (370.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH1'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.430d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH1'
      vnlocal(isame,2) = 'CT1'
c     --- NH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH1'
      vnlocal(isame,2) = 'CT2'
c     --- NH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH1'
      vnlocal(isame,2) = 'CT3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.997d0
      vclocal(1) = (440.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH1'
      vnlocal(isame,2) = 'H'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9800d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH1 - HC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH1'
      vnlocal(isame,2) = 'HC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3600d0
      vclocal(1) = (430.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH2'
      vnlocal(isame,2) = 'CC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4550d0
      vclocal(1) = (240.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH2'
      vnlocal(isame,2) = 'CT2'
c     --- NH2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH2'
      vnlocal(isame,2) = 'CT3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (480.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH2 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH2'
      vnlocal(isame,2) = 'H'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (460.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH2 - HC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH2'
      vnlocal(isame,2) = 'HC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4800d0
      vclocal(1) = (200.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH3 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH3'
      vnlocal(isame,2) = 'CT1'
c     --- NH3 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH3'
      vnlocal(isame,2) = 'CT2'
c     --- NH3 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH3'
      vnlocal(isame,2) = 'CT3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0400d0
      vclocal(1) = (403.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NH3 - HC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NH3'
      vnlocal(isame,2) = 'HC'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4850d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NP - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NP'
      vnlocal(isame,2) = 'CP1'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5020d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NP - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NP'
      vnlocal(isame,2) = 'CP3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0060d0
      vclocal(1) = (460.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NP - HC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NP'
      vnlocal(isame,2) = 'HC'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3757d0
      vclocal(1) = (377.2d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NPH - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NPH'
      vnlocal(isame,2) = 'CPA'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3800d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NR1'
      vnlocal(isame,2) = 'CPH1'
c     --- NR2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NR2'
      vnlocal(isame,2) = 'CPH1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3600d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NR1'
      vnlocal(isame,2) = 'CPH2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3200d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR2 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NR2'
      vnlocal(isame,2) = 'CPH2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (466.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NR1'
      vnlocal(isame,2) = 'H'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3700d0
      vclocal(1) = (380.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR3 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NR3'
      vnlocal(isame,2) = 'CPH1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3200d0
      vclocal(1) = (380.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR3 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NR3'
      vnlocal(isame,2) = 'CPH2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (453.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR3 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NR3'
      vnlocal(isame,2) = 'H'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3700d0
      vclocal(1) = (270.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NY'
      vnlocal(isame,2) = 'CA'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3750d0
      vclocal(1) = (270.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NY'
      vnlocal(isame,2) = 'CPT'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9760d0
      vclocal(1) = (465.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NY - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NY'
      vnlocal(isame,2) = 'H'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.2300d0
      vclocal(1) = (620.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- O - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'O'
      vnlocal(isame,2) = 'C'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.2300d0
      vclocal(1) = (650.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- O - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'O'
      vnlocal(isame,2) = 'CC'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3300d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OC - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OC'
      vnlocal(isame,2) = 'CT2'
c     --- OC - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OC'
      vnlocal(isame,2) = 'CT3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4110d0
      vclocal(1) = (334.3d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OH1 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OH1'
      vnlocal(isame,2) = 'CA'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.1280d0
      vclocal(1) = (1115.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OM - CM [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OM'
      vnlocal(isame,2) = 'CM'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.2300d0
      vclocal(1) = (600.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OM - OM [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OM'
      vnlocal(isame,2) = 'OM'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.8180d0
      vclocal(1) = (198.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- S - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'S'
      vnlocal(isame,2) = 'CT2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.8160d0
      vclocal(1) = (240.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- S - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'S'
      vnlocal(isame,2) = 'CT3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3250d0
      vclocal(1) = (275.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- S - HS [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'S'
      vnlocal(isame,2) = 'HS'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.8160d0
      vclocal(1) = (214.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- SM - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'SM'
      vnlocal(isame,2) = 'CT2'
c     --- SM - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'SM'
      vnlocal(isame,2) = 'CT3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.0290d0
      vclocal(1) = (173.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- SM - SM [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'SM'
      vnlocal(isame,2) = 'SM'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.8360d0
      vclocal(1) = (205.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- SS - CS [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'SS'
      vnlocal(isame,2) = 'CS'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     *** Nucleic acid bond terms (if not listed above) ***
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.090d0
      vclocal(1) = (374.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3C - HN6 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3C'
      vnlocal(isame,2) = 'HN6'
c     --- CN3C - HN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3C'
      vnlocal(isame,2) = 'HN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.090d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - HN6 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3'
      vnlocal(isame,2) = 'HN6'
c     --- CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3'
      vnlocal(isame,2) = 'HN3'
c     --- CN3T - HN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3T'
      vnlocal(isame,2) = 'HN3'
c     --- CN3 - HN3B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3'
      vnlocal(isame,2) = 'HN3B'
c     --- CN3A - HN3B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3A'
      vnlocal(isame,2) = 'HN3B'
c     --- CN3B - HN3B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3B'
      vnlocal(isame,2) = 'HN3B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.409d0
      vclocal(1) = (302.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - CN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'CN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.403d0
      vclocal(1) = (302.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - CN3T [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'CN3T'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.480d0
      vclocal(1) = (302.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1A - CN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1A'
      vnlocal(isame,2) = 'CN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.360d0
      vclocal(1) = (302.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - CN5G [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'CN5G'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.450d0
      vclocal(1) = (200.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - CN7B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'CN7B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.518d0
      vclocal(1) = (200.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - CN7C [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'CN7C'
c     --- CN7B - CN8 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'CN8'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.456d0
      vclocal(1) = (220.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'NN2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.458d0
      vclocal(1) = (220.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'NN2B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.462d0
      vclocal(1) = (220.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'NR1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.460d0
      vclocal(1) = (222.5d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN7B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'CN7B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.490d0
      vclocal(1) = (222.5d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - CN8 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3'
      vnlocal(isame,2) = 'CN8'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.512d0
      vclocal(1) = (222.5d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'CN8B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.516d0
      vclocal(1) = (222.5d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN7C [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'CN7C'
c     --- CN7 - CN8 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'CN8'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.529d0
      vclocal(1) = (222.5d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'CN7'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.478d0
      vclocal(1) = (230.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3T - CN9 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3T'
      vnlocal(isame,2) = 'CN9'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.480d0
      vclocal(1) = (230.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3D - CN9 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3D'
      vnlocal(isame,2) = 'CN9'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.580d0
      vclocal(1) = (237.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- ON4 - P [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'ON4'
      vnlocal(isame,2) = 'P'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.446d0
      vclocal(1) = (240.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'ON6'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.480d0
      vclocal(1) = (240.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - ON6B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'ON6B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4050d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - CN3D [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN2'
      vnlocal(isame,2) = 'CN3D'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.420d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - ON6B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'ON6B'
c     --- CN7B - ON6 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'ON6'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.375d0
      vclocal(1) = (300.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5 - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN5'
      vnlocal(isame,2) = 'NN2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.378d0
      vclocal(1) = (300.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN4 - NN2B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN4'
      vnlocal(isame,2) = 'NN2B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.680d0
      vclocal(1) = (300.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - P2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'ON2'
      vnlocal(isame,2) = 'P2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.343d0
      vclocal(1) = (302.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3'
      vnlocal(isame,2) = 'NN2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.348d0
      vclocal(1) = (302.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1T - NN2B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1T'
      vnlocal(isame,2) = 'NN2B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.375d0
      vclocal(1) = (302.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5 - NN2B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN5'
      vnlocal(isame,2) = 'NN2B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.355d0
      vclocal(1) = (310.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5 - NN4 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN5'
      vnlocal(isame,2) = 'NN4'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.361d0
      vclocal(1) = (310.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5 - CN5 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN5'
      vnlocal(isame,2) = 'CN5'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (310.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5G - NN4 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN5G'
      vnlocal(isame,2) = 'NN4'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.433d0
      vclocal(1) = (310.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7'
      vnlocal(isame,2) = 'ON2'
c     --- CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'ON2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.326d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - NN3G [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN2'
      vnlocal(isame,2) = 'NN3G'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.343d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - NN2B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3'
      vnlocal(isame,2) = 'NN2B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.350d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5 - CN5G [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN5'
      vnlocal(isame,2) = 'CN5G'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.374d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN4 - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN4'
      vnlocal(isame,2) = 'NN2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.406d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN2'
      vnlocal(isame,2) = 'CN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.440d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8B - ON2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8B'
      vnlocal(isame,2) = 'ON2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.383d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1T - NN2U [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1T'
      vnlocal(isame,2) = 'NN2U'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.389d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - NN2U [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'NN2U'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.396d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - NN2G [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'NN2G'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.440d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8 - ON2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8'
      vnlocal(isame,2) = 'ON2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.312d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5 - NN3A [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN5'
      vnlocal(isame,2) = 'NN3A'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.315d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5 - NN3G [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN5'
      vnlocal(isame,2) = 'NN3G'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.332d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5 - NN3I [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN5'
      vnlocal(isame,2) = 'NN3I'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - NN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'NN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.358d0
      vclocal(1) = (360.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - CN5 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN2'
      vnlocal(isame,2) = 'CN5'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.366d0
      vclocal(1) = (360.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - NN1 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN2'
      vnlocal(isame,2) = 'NN1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.090d0
      vclocal(1) = (380.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN4 - HN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN4'
      vnlocal(isame,2) = 'HN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.367d0
      vclocal(1) = (380.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'NN2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.305d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN4 - NN4 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN4'
      vnlocal(isame,2) = 'NN4'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.342d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - NN3A [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN2'
      vnlocal(isame,2) = 'NN3A'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.392d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - NN2G [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN2'
      vnlocal(isame,2) = 'NN2G'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.456d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN9 - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN9'
      vnlocal(isame,2) = 'NN2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.458d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8 - NN2B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8'
      vnlocal(isame,2) = 'NN2B'
c     --- CN9 - NN2B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN9'
      vnlocal(isame,2) = 'NN2B'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.460d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0 
c     --- CN8 - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN8'
      vnlocal(isame,2) = 'NN2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.315d0
      vclocal(1) = (420.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3B - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3B'
      vnlocal(isame,2) = 'NN2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.320d0
      vclocal(1) = (420.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3C - CN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3C'
      vnlocal(isame,2) = 'CN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.322d0
      vclocal(1) = (420.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN4 - NN3A [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN4'
      vnlocal(isame,2) = 'NN3A'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.350d0
      vclocal(1) = (420.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3B - CN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3B'
      vnlocal(isame,2) = 'CN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.355d0
      vclocal(1) = (420.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3C - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3C'
      vnlocal(isame,2) = 'NN2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.400d0
      vclocal(1) = (428.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - ON5 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'ON5'
c     --- CN7C - ON5 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN7C'
      vnlocal(isame,2) = 'ON5'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.295d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN4 - NN3I [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN4'
      vnlocal(isame,2) = 'NN3I'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.343d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - NN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN2'
      vnlocal(isame,2) = 'NN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.360d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3A - CN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3A'
      vnlocal(isame,2) = 'CN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.010d0
      vclocal(1) = (460.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NN5 - HN1 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NN5'
      vnlocal(isame,2) = 'HN1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.010d0
      vclocal(1) = (471.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HN2 - NN2G [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HN2'
      vnlocal(isame,2) = 'NN2G'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.010d0
      vclocal(1) = (474.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HN2 - NN2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HN2'
      vnlocal(isame,2) = 'NN2'
c     --- HN2 - NN2B [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HN2'
      vnlocal(isame,2) = 'NN2B'
c     --- HN2 - NN2U [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HN2'
      vnlocal(isame,2) = 'NN2U'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.530d0
      vclocal(1) = (480.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- ON3 - P2 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'ON3'
      vnlocal(isame,2) = 'P2'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.000d0
      vclocal(1) = (488.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- HN1 - NN1 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'HN1'
      vnlocal(isame,2) = 'NN1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.326d0
      vclocal(1) = (500.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - CN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3'
      vnlocal(isame,2) = 'CN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.320d0
      vclocal(1) = (560.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - CN3T [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3'
      vnlocal(isame,2) = 'CN3T'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (560.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3D - CN3 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN3D'
      vnlocal(isame,2) = 'CN3'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.360d0
      vclocal(1) = (560.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1A - NN1 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1A'
      vnlocal(isame,2) = 'NN1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.245d0
      vclocal(1) = (620.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - ON1C [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'ON1C'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.234d0
      vclocal(1) = (660.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - ON1 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'ON1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.230d0
      vclocal(1) = (860.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1A - ON1 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1A'
      vnlocal(isame,2) = 'ON1'
c     --- CN1T - ON1 [par_all27_na.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1T'
      vnlocal(isame,2) = 'ON1'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.9000d0
      vclocal(1) = (258.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- FE - CM [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'FE'
      vnlocal(isame,2) = 'CM'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 3.3814d0
      vclocal(1) = (0.000d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- FE - CPM [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'FE'
      vnlocal(isame,2) = 'CPM'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.9580d0
      vclocal(1) = (270.200d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NPH - FE [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NPH'
      vnlocal(isame,2) = 'FE'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.200d0
      vclocal(1) = (65.000d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR2 - FE [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'NR2'
      vnlocal(isame,2) = 'FE'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.800d0
      vclocal(1) = (250.000d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- OM - FE [par_all27_prot_lipid.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'OM'
      vnlocal(isame,2) = 'FE'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return
      return
      end

      subroutine twh_ffcharmm27_angle(lfailure,ffname,nangletypes)
c     ******************************************************************
c     * sets up the angle terms for Charmm27                           *
c     *                                                                *
c     * split from original ffcharmm27 03-05-2006 by M.G. Martin       *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nangletypes
c     --- local variables
      integer type,isame,nbslocal
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)

c     --- initialize type count
      type = 0
c     ******************************************************************
c     * charmm27 angles with urey-bradley                              *
c     ******************************************************************
      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (52.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CTL2 - CTL2 - CL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CL'
C     --- CTL3 - CTL2 - CL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL3'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CL'
C     --- CT1 - CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'C'
C     --- CT1 - CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CC'
C     --- CT1 - CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CC'
C     --- CT1 - CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CD'
C     --- CT2 - CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'C'
C     --- CT2 - CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CC'
C     --- CT2 - CT1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CD'
C     --- CT2 - CT2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'C'
C     --- CT2 - CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CC'
C     --- CT2 - CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CD'
C     --- CT3 - CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'C'
C     --- CT3 - CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CC'
C     --- NH2 - CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.00d0
      bclocal(1) = (53.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CTL1 - CTL1 - CTL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL1'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL1'
C     --- CT1 - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT1'
C     --- CT2 - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 113.50d0
      bclocal(1) = (58.350d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (11.16d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CTL1 - CTL1 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL1'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL2'
C     --- CTL1 - CTL2 - CTL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL1'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL1'
C     --- CTL1 - CTL2 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL1'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL2'
C     --- CTL1 - CTL2 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL1'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL3'
C     --- CTL2 - CTL1 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL2'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL2'
C     --- CTL2 - CTL1 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL2'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL3'
C     --- CT1 - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT1'
C     --- CT2 - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT1'
C     --- CT2 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
C     --- CT2 - CT3 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CT1'
C     --- CT3 - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 113.60d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (11.16d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CTL2 - CTL2 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 113.60d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (11.16d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7B - CN7 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN8'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.60d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2651d0
      bclocal(3) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CL - OSL - CTL1 methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CL'
      bnlocal(isame,2) = 'OSL'
      bnlocal(isame,3) = 'CTL1'
C     --- CL - OSL - CTL2 methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CL'
      bnlocal(isame,2) = 'OSL'
      bnlocal(isame,3) = 'CTL2'
C     --- CL - OSL - CTL3 methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CL'
      bnlocal(isame,2) = 'OSL'
      bnlocal(isame,3) = 'CTL3'
C     --- CT2 - OS - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'CD'
C     --- CT3 - OS - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'OS'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.50d0
      bclocal(1) = (53.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CTL1 - CTL1 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL1'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL3'
C     --- CT3 - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 115.00d0
      bclocal(1) = (58.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CTL2 - CTL2 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL3'
C     --- CT3 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.00d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HOL - OHL - CL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HOL'
      bnlocal(isame,2) = 'OHL'
      bnlocal(isame,3) = 'CL'
C     --- H - OH1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'OH1'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.00d0
      bclocal(1) = (57.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HOL - OHL - CTL1
      isame = isame + 1
      bnlocal(isame,1) = 'HOL'
      bnlocal(isame,2) = 'OHL'
      bnlocal(isame,3) = 'CTL1'
C     --- HOL - OHL - CTL2
      isame = isame + 1
      bnlocal(isame,1) = 'HOL'
      bnlocal(isame,2) = 'OHL'
      bnlocal(isame,3) = 'CTL2'
C     --- HOL - OHL - CTL3
      isame = isame + 1
      bnlocal(isame,1) = 'HOL'
      bnlocal(isame,2) = 'OHL'
      bnlocal(isame,3) = 'CTL3'
C     --- H - OH1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'OH1'
      bnlocal(isame,3) = 'CT1'
C     --- H - OH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'OH1'
      bnlocal(isame,3) = 'CT2'
C     --- H - OH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'OH1'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.00d0
      bclocal(1) = (57.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN5 - ON5 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN5'
      bnlocal(isame,2) = 'ON5'
      bnlocal(isame,3) = 'CN8'
C     --- HN5 - ON5 - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN5'
      bnlocal(isame,2) = 'ON5'
      bnlocal(isame,3) = 'CN9'
C     --- HN5 - ON5 - CN8B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN5'
      bnlocal(isame,2) = 'ON5'
      bnlocal(isame,3) = 'CN8B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (33.430d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HL - CTL2 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL2'
C     --- HL - CTL2 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL3'
C     --- HAL3 - CTL3 - CTL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL3'
      bnlocal(isame,2) = 'CTL3'
      bnlocal(isame,3) = 'CTL1'
C     --- HA - CP2 - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CP2'
      bnlocal(isame,3) = 'CP1'
C     --- HA - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT1'
C     --- HA - CT3 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CT1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (33.430d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN9'
      bnlocal(isame,3) = 'HN9'
C     --- CN3T - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3T'
      bnlocal(isame,2) = 'CN9'
      bnlocal(isame,3) = 'HN9'
C     --- CN3D - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3D'
      bnlocal(isame,2) = 'CN9'
      bnlocal(isame,3) = 'HN9'
C     --- CN7B - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'HN8'
C     --- CN8 - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'HN7'
C     --- CN7B - CN7C - HN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7C'
      bnlocal(isame,3) = 'HN7'
C     --- CN7C - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7C'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'HN7'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (26.50d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HAL2 - CTL2 - CTL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL1'
C     --- HAL2 - CTL2 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL2'
C     --- HA - CP2 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CP2'
      bnlocal(isame,3) = 'CP2'
C     --- HA - CP2 - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CP2'
      bnlocal(isame,3) = 'CP3'
C     --- HA - CP3 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CP3'
      bnlocal(isame,3) = 'CP2'
C     --- HA - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.00d0
      bclocal(1) = (35.50d0)*(CON_KCALTOK)
      bclocal(2) = 1.8020d0
      bclocal(3) = (5.40d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HAL2 - CTL2 - HAL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'HAL2'
C     --- HA - CP2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CP2'
      bnlocal(isame,3) = 'HA'
C     --- HA - CP3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CP3'
      bnlocal(isame,3) = 'HA'
C     --- HA - CT1 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'HA'
C     --- HA - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.00d0
      bclocal(1) = (35.50d0)*(CON_KCALTOK)
      bclocal(2) = 1.8020d0
      bclocal(3) = (5.40d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN8 - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'HN8'
C     --- HN7 - CN7 - HN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'HN7'
C     --- HN7 - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'HN7'
C     --- HN7 - CN7C - HN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7C'
      bnlocal(isame,3) = 'HN7'
C     --- HN8 - CN8B - HN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8B'
      bnlocal(isame,3) = 'HN8'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (34.60d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HAL2 - CTL2 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL3'
C     --- HAL3 - CTL3 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL3'
      bnlocal(isame,2) = 'CTL3'
      bnlocal(isame,3) = 'CTL2'
C     --- HA - CS - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CS'
      bnlocal(isame,3) = 'CT3'
C     --- HA - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT3'
C     --- HA - CT3 - CS [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CS'
C     --- HA - CT3 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (34.60d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN9 - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN9'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'HN8'
C     --- CN8 - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8'
      bnlocal(isame,2) = 'CN9'
      bnlocal(isame,3) = 'HN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (33.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1630d0
      bclocal(3) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HAL2 - CTL2 - CL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CL'
C     --- HAL3 - CTL3 - CL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL3'
      bnlocal(isame,2) = 'CTL3'
      bnlocal(isame,3) = 'CL'
C     --- HA - CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'C'
C     --- HA - CT1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CD'
C     --- HA - CT2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'C'
C     --- HA - CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CC'
C     --- HA - CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CD'
C     --- HA - CT3 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'C'
C     --- HA - CT3 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CC'
C     --- HA - CT3 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame 
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (34.50d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HAL1 - CTL1 - CTL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL1'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL1'
C     --- HAL1 - CTL1 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL1'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL2'
C     --- HAL1 - CTL1 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL1'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL3'
C     --- HA - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT1'
C     --- HA - CT1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT2'
C     --- HA - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (34.50d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN7 - CN7 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN8'
C     --- HN7 - CN7 - CN7C [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7C'
C     --- HN7 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN8B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (37.50d0)*(CON_KCALTOK)
      bclocal(2) = 2.1790d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HAL3 - CTL3 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL3'
      bnlocal(isame,2) = 'CTL3'
      bnlocal(isame,3) = 'CTL3'
C     --- HA - CT3 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.40d0
      bclocal(1) = (35.50d0)*(CON_KCALTOK)
      bclocal(2) = 1.8020d0
      bclocal(3) = (5.40d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HAL3 - CTL3 - HAL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL3'
      bnlocal(isame,2) = 'CTL3'
      bnlocal(isame,3) = 'HAL3'
C     --- HA - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
C     --- HN9 - CN9 - HN9 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN9'
      bnlocal(isame,2) = 'CN9'
      bnlocal(isame,3) = 'HN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (33.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.0560d0
      bclocal(3) = (4.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HCL - NH3L - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HCL'
      bnlocal(isame,2) = 'NH3L'
      bnlocal(isame,3) = 'CTL2'
C     --- HC - NP - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NP'
      bnlocal(isame,3) = 'CP1'
C     --- HC - NP - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NP'
      bnlocal(isame,3) = 'CP3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.52d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HT - OT - HT [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HT'
      bnlocal(isame,2) = 'OT'
      bnlocal(isame,3) = 'HT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.00d0
      bclocal(1) = (67.70d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH3L - CTL2 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3L'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL2'
C     --- NH3 - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT1'
C     --- NH3 - CT1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT2'
C     --- NH3 - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT3'
C     --- NH3 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.00d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.4420d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OBL - CL - CTL3 methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OBL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'CTL3'
C     --- OBL - CL - CTL2 methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OBL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'CTL2'
C     --- OBL - CL - CTL1 methyl acetate [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OBL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'CTL1'
C     --- OB - CD - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OB'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CP1'
C     --- OB - CD - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OB'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CT1'
C     --- OB - CD - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OB'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CT2'
C     --- OB - CD - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OB'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 118.00d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.3880d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OCL - CL - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OCL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'CTL2'
C     --- OCL - CL - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OCL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'CTL3'
C     --- OC - CC - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OC'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CP1'
C     --- OC - CC - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OC'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT1'
C     --- OC - CC - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OC'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT2'
C     --- OC - CC - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OC'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 124.00d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2250d0
      bclocal(3) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OCL - CL - OCL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OCL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'OCL'
C     --- OC - CC - OC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OC'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'OC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.50d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OHL - CL - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'CTL3'
C     --- OHL - CL - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'CTL2'
C     --- OH1 - CD - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CT2'
C     --- OH1 - CD - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 123.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2620d0
      bclocal(3) = (210.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OHL - CL - OBL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'OBL'
C     --- OH1 - CD - OB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'OB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.10d0
      bclocal(1) = (75.70d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OSL - CTL1 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL2'
C     --- OSL - CTL1 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL3'
C     --- OSL - CTL2 - CTL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL1'
C     --- OSL - CTL2 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL2'
C     --- OSL - CTL2 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL3'
C     --- OHL - CTL1 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'CTL2'
C     --- OHL - CTL2 - CTL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL1'
C     --- OHL - CTL2 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL2'
C     --- OHL - CTL2 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL3'
C     --- OH1 - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT1'
C     --- OH1 - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT3'
C     --- OH1 - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT1'
C     --- OH1 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
C     --- OH1 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.10d0
      bclocal(1) = (75.70d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7 - CN8 - ON5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'ON5'
C     --- CN7 - CN8B - ON5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'CN8B'
      bnlocal(isame,3) = 'ON5'
C     --- ON5 - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON5'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7'
C     --- ON5 - CN7 - CN7C [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON5'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.89d0
      bclocal(1) = (45.90d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OHL - CTL1 - HAL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'HAL1'
C     --- OHL - CTL2 - HAL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'HAL2'
C     --- OHL - CTL3 - HAL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OHL'
      bnlocal(isame,2) = 'CTL3'
      bnlocal(isame,3) = 'HAL3'
C     --- OH1 - CT1 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'HA'
C     --- OH1 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA'
C     --- OH1 - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.89d0
      bclocal(1) = (45.90d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN8 - CN8 - ON5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'ON5'
C     --- HN8 - CN8B - ON5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8B'
      bnlocal(isame,3) = 'ON5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.00d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.3260d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OSL - CL - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'CTL3'
C     --- OSL - CL - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'CTL2'
C     --- OS - CD - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CP1'
C     --- OS - CD - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CT1'
C     --- OS - CD - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CT2'
C     --- OS - CD - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.90d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2576d0
      bclocal(3) = (160.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OSL - CL - OBL acetic acid [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'CL'
      bnlocal(isame,3) = 'OBL'
C     --- OS - CD - OB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'CD'
      bnlocal(isame,3) = 'OB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HAL1 - CTL1 - OSL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL1'
      bnlocal(isame,2) = 'CTL1'
      bnlocal(isame,3) = 'OSL'
C     --- HAL2 - CTL2 - OSL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'OSL'
C     --- HAL3 - CTL3 - OSL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL3'
      bnlocal(isame,2) = 'CTL3'
      bnlocal(isame,3) = 'OSL'
C     --- OS - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA'
C     --- OS - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OS'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN9 - CN9 - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN9'
      bnlocal(isame,2) = 'CN9'
      bnlocal(isame,3) = 'ON2'
C     --- HN8 - CN8B - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8B'
      bnlocal(isame,3) = 'ON2'
C     --- HN7 - CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'ON2'
C     --- HN7 - CN7B - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'ON2'
C     --- HN7 - CN7 - ON5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'ON5'
C     --- HN7 - CN7B - ON5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'ON5'
C     --- HN7 - CN7C - ON5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7C'
      bnlocal(isame,3) = 'ON5'
C     --- HN8 - CN8 - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'ON2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.0d0
      bclocal(1) = (20.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.33d0
      bclocal(3) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL2 - OSL - PL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL2'
      bnlocal(isame,2) = 'OSL'
      bnlocal(isame,3) = 'PL'
c     --- CTL3 - OSL - PL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL3'
      bnlocal(isame,2) = 'OSL'
      bnlocal(isame,3) = 'PL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.0d0
      bclocal(1) = (20.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.33d0
      bclocal(3) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CN9 - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN9'
      bnlocal(isame,2) = 'ON2'
      bnlocal(isame,3) = 'P'
c     --- CN9 - ON2 - P2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN9'
      bnlocal(isame,2) = 'ON2'
      bnlocal(isame,3) = 'P2'
c     --- CN8B - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8B'
      bnlocal(isame,2) = 'ON2'
      bnlocal(isame,3) = 'P'
c     --- CN7B - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'ON2'
      bnlocal(isame,3) = 'P'
c     --- CN7 - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'ON2'
      bnlocal(isame,3) = 'P'
c     --- CN8 - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8'
      bnlocal(isame,2) = 'ON2'
      bnlocal(isame,3) = 'P'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 115.0d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.30d0
      bclocal(3) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HOL - OHL - PL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HOL'
      bnlocal(isame,2) = 'OHL'
      bnlocal(isame,3) = 'PL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.3d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- OSL - PL - OSL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'PL'
      bnlocal(isame,3) = 'OSL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.3d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - P - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON2'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'ON2'
c     --- ON2 - P2 - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON2'
      bnlocal(isame,2) = 'P2'
      bnlocal(isame,3) = 'ON2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.6d0
      bclocal(1) = (98.9d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- OSL - PL - O2L [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'PL'
      bnlocal(isame,3) = 'O2L'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.6d0
      bclocal(1) = (98.9d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - P - ON3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON2'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'ON3'
c     --- ON5 - P - ON3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON5'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'ON3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.0d0
      bclocal(1) = (48.1d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- OSL - PL - OHL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OSL'
      bnlocal(isame,2) = 'PL'
      bnlocal(isame,3) = 'OHL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.0d0
      bclocal(1) = (48.1d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - P - ON4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON2'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'ON4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (120.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- O2L - PL - O2L [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O2L'
      bnlocal(isame,2) = 'PL'
      bnlocal(isame,3) = 'O2L'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (120.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- ON3 - P - ON3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON3'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'ON3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.23d0
      bclocal(1) = (98.9d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- O2L - PL - OHL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O2L'
      bnlocal(isame,2) = 'PL'
      bnlocal(isame,3) = 'OHL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.23d0
      bclocal(1) = (98.9d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- ON3 - P - ON4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON3'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'ON4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.5d0
      bclocal(1) = (40.0d0)*(CON_KCALTOK)
      bclocal(2) = 2.13d0
      bclocal(3) = (27.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- NTL - CTL2 - HL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NTL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'HL'
c     --- NTL - CTL5 - HL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NTL'
      bnlocal(isame,2) = 'CTL5'
      bnlocal(isame,3) = 'HL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (24.0d0)*(CON_KCALTOK)
      bclocal(2) = 1.767d0
      bclocal(3) = (28.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HL - CTL2 - HL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'HL'
c     --- HL - CTL5 - HL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HL'
      bnlocal(isame,2) = 'CTL5'
      bnlocal(isame,3) = 'HL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.5d0
      bclocal(1) = (60.0d0)*(CON_KCALTOK)
      bclocal(2) = 2.466d0
      bclocal(3) = (26.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL5 - NTL - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL5'
      bnlocal(isame,2) = 'NTL'
      bnlocal(isame,3) = 'CTL2'
c     --- CTL5 - NTL - CTL5 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL5'
      bnlocal(isame,2) = 'NTL'
      bnlocal(isame,3) = 'CTL5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.00d0
      bclocal(1) = (67.7d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- NTL - CTL2 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NTL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL2'
c     --- NTL - CTL2 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NTL'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (41.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HCL - NH3L - HCL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HCL'
      bnlocal(isame,2) = 'NH3L'
      bnlocal(isame,3) = 'HCL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.50d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      bclocal(2) = 2.0836d0
      bclocal(3) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- NH3L - CTL2 - HAL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3L'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'HAL2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.47d0
      bclocal(1) = (130.0d0)*(CON_KCALTOK)
      bclocal(2) = 2.45d0
      bclocal(3) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- O2L - SL - O2L [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O2L'
      bnlocal(isame,2) = 'SL'
      bnlocal(isame,3) = 'O2L'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 98.0d0
      bclocal(1) = (85.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- O2L - SL - OSL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O2L'
      bnlocal(isame,2) = 'SL'
      bnlocal(isame,3) = 'OSL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.0d0
      bclocal(1) = (15.0d0)*(CON_KCALTOK)
      bclocal(2) = 1.90d0
      bclocal(3) = (27.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL2 - OSL - SL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL2'
      bnlocal(isame,2) = 'OSL'
      bnlocal(isame,3) = 'SL'
c     --- CTL3 - OSL - SL [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CTL3'
      bnlocal(isame,2) = 'OSL'
      bnlocal(isame,3) = 'SL'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.50d0
      bclocal(1) = (48.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL1 - CEL1 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CEL1'
      bnlocal(isame,2) = 'CEL1'
      bnlocal(isame,3) = 'CTL2'
c     --- CEL1 - CEL1 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CEL1'
      bnlocal(isame,2) = 'CEL1'
      bnlocal(isame,3) = 'CTL3'
c     --- CE1 - CE1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CE1'
      bnlocal(isame,2) = 'CE1'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.00d0
      bclocal(1) = (48.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL2 - CEL1 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CEL2'
      bnlocal(isame,2) = 'CEL1'
      bnlocal(isame,3) = 'CTL2'
c     --- CE2 - CE1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CE2'
      bnlocal(isame,2) = 'CE1'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.20d0
      bclocal(1) = (47.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL2 - CEL1 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CEL2'
      bnlocal(isame,2) = 'CEL1'
      bnlocal(isame,3) = 'CTL3'
c     --- CE2 - CE1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CE2'
      bnlocal(isame,2) = 'CE1'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.50d0
      bclocal(1) = (52.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL1 - CEL1 - CEL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HEL1'
      bnlocal(isame,2) = 'CEL1'
      bnlocal(isame,3) = 'CEL1'
c     --- HE1 - CE1 - CE1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HE1'
      bnlocal(isame,2) = 'CE1'
      bnlocal(isame,3) = 'CE1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (42.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL1 - CEL1 - CEL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HEL1'
      bnlocal(isame,2) = 'CEL1'
      bnlocal(isame,3) = 'CEL2'
c     --- HE1 - CE1 - CE2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HE1'
      bnlocal(isame,2) = 'CE1'
      bnlocal(isame,3) = 'CE2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.00d0
      bclocal(1) = (40.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL1 - CEL1 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HEL1'
      bnlocal(isame,2) = 'CEL1'
      bnlocal(isame,3) = 'CTL2'
c     --- HE1 - CE1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HE1'
      bnlocal(isame,2) = 'CE1'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (22.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL1 - CEL1 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HEL1'
      bnlocal(isame,2) = 'CEL1'
      bnlocal(isame,3) = 'CTL3'
c     --- HE1 - CE1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HE1'
      bnlocal(isame,2) = 'CE1'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL2 - CEL2 - CEL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HEL2'
      bnlocal(isame,2) = 'CEL2'
      bnlocal(isame,3) = 'CEL1'
c     --- HE2 - CE2 - CE1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HE2'
      bnlocal(isame,2) = 'CE2'
      bnlocal(isame,3) = 'CE1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (55.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL2 - CEL2 - CEL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HEL2'
      bnlocal(isame,2) = 'CEL2'
      bnlocal(isame,3) = 'CEL2'
c     --- HE2 - CE2 - CE2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HE2'
      bnlocal(isame,2) = 'CE2'
      bnlocal(isame,3) = 'CE2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.00d0
      bclocal(1) = (19.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL2 - CEL2 - HEL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HEL2'
      bnlocal(isame,2) = 'CEL2'
      bnlocal(isame,3) = 'HEL2'
c     --- HE2 - CE2 - HE2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HE2'
      bnlocal(isame,2) = 'CE2'
      bnlocal(isame,3) = 'HE2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.20d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL1 - CTL2 - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CEL1'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL2'
c     --- CEL1 - CTL2 - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CEL1'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CTL3'
c     --- CE1 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CE1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.50d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HAL2 - CTL2 - CEL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL2'
      bnlocal(isame,2) = 'CTL2'
      bnlocal(isame,3) = 'CEL1'
c     --- HA - CT2 - CE1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CE1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.50d0
      bclocal(1) = (42.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- HAL3 - CTL3 - CEL1 [par_all27_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HAL3'
      bnlocal(isame,2) = 'CTL3'
      bnlocal(isame,3) = 'CEL1'
c     --- HA - CT3 - CE1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CE1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.4162d0
      bclocal(3) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CA - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CP1 - N - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CP1'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'C'
c     --- CP3 - N - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CP3'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly) 
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.30d0
      bclocal(1) = (52.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CP2 - CP1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CP2'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'C'
c     --- CP2 - CP1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CP2'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.30d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CP2 - CP1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CP2'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CP2 - CP2 - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CP2'
      bnlocal(isame,2) = 'CP2'
      bnlocal(isame,3) = 'CP1'
c     --- CP3 - CP2 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CP3'
      bnlocal(isame,2) = 'CP2'
      bnlocal(isame,3) = 'CP2'
c     --- NP - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NP'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CP2'
c     --- NP - CP3 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NP'
      bnlocal(isame,2) = 'CP3'
      bnlocal(isame,3) = 'CP2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.20d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CP3 - N - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CP3'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'CP1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.00d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CP3 - NP - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CP3'
      bnlocal(isame,2) = 'NP'
      bnlocal(isame,3) = 'CP1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.74d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPA - CPB - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPA'
      bnlocal(isame,2) = 'CPB'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.12d0
      bclocal(1) = (94.20d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPA - CPM - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPA'
      bnlocal(isame,2) = 'CPM'
      bnlocal(isame,3) = 'CPA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.90d0
      bclocal(1) = (139.30d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPA - NPH - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPA'
      bnlocal(isame,2) = 'NPH'
      bnlocal(isame,3) = 'CPA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPB - C - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPB'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.75d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPB - CPB - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPB'
      bnlocal(isame,2) = 'CPB'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.51d0
      bclocal(1) = (30.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPB - CPB - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPB'
      bnlocal(isame,2) = 'CPB'
      bnlocal(isame,3) = 'CPA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.50d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH2 - NR1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPH2'
      bnlocal(isame,2) = 'NR1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.00d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH2 - NR2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPH2'
      bnlocal(isame,2) = 'NR2'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (145.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH2 - NR3 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPH2'
      bnlocal(isame,2) = 'NR3'
      bnlocal(isame,3) = 'CPH1'
c     --- NR3 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR3'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CPH1'
c     --- NR3 - CPH2 - NR3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR3'
      bnlocal(isame,2) = 'CPH2'
      bnlocal(isame,3) = 'NR3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.07d0
      bclocal(1) = (61.60d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPM - CPA - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPM'
      bnlocal(isame,2) = 'CPA'
      bnlocal(isame,3) = 'CPB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPT - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPT'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPT - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPT'
      bnlocal(isame,2) = 'CPT'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - CN3C - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3C'
      bnlocal(isame,3) = 'NN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.40d0
      bclocal(1) = (120.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.261d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPT - CY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPT'
      bnlocal(isame,2) = 'CY'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (110.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CPT - NY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPT'
      bnlocal(isame,2) = 'NY'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (110.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- ON6 - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2'
c     --- ON6 - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2B'
c     --- CN7 - ON6 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'ON6'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.50d0
      bclocal(1) = (51.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT1 - CT2 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CA'
c     --- CT3 - CT2 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.00d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT1 - CT2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CPH1'
c     --- CT3 - CT2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT1 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT1'
      bnlocal(isame,2) = 'NH1'
      bnlocal(isame,3) = 'C'
c     --- CT2 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'NH1'
      bnlocal(isame,3) = 'C'
c     --- CT3 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'NH1'
      bnlocal(isame,3) = 'C'
c     --- H - NH2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NH2'
      bnlocal(isame,3) = 'CC'
c     --- HA - C - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CPB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.30d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     --- CT3 - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.74d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CPB - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CPB'
      bnlocal(isame,3) = 'CPA'
c     --- CT3 - CPB - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CPB'
      bnlocal(isame,3) = 'CPA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.75d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CPB - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CPB'
      bnlocal(isame,3) = 'CPB'
c     --- CT3 - CPB - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CPB'
      bnlocal(isame,3) = 'CPB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.00d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CPH1'
c     --- CT3 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.00d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CT2 - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CPB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.40d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CY'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.00d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'CY'
      bnlocal(isame,3) = 'CPT'
c     --- NR1 - CPH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR1'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CT2'
c     --- NR1 - CPH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR1'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (62.30d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - NC2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT2'
      bnlocal(isame,2) = 'NC2'
      bnlocal(isame,3) = 'C'
c     --- CT3 - NC2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'NC2'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 114.00d0
      bclocal(1) = (53.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.5610d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CT3 - CT1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT2'
C     --- CT3 - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT3'
C     --- CT3 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 95.00d0
      bclocal(1) = (34.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CT3 - S - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.60d0
      bclocal(1) = (160.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CY - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CY'
      bnlocal(isame,2) = 'CPT'
      bnlocal(isame,3) = 'CA'
C     --- NY - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NY'
      bnlocal(isame,2) = 'CPT'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.40d0
      bclocal(1) = (110.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CY - CPT - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CY'
      bnlocal(isame,2) = 'CPT'
      bnlocal(isame,3) = 'CPT'
C     --- NY - CPT - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NY'
      bnlocal(isame,2) = 'CPT'
      bnlocal(isame,3) = 'CPT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.00d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CY - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CY'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT1'
C     --- CY - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CY'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.00d0
      bclocal(1) = (34.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NH1'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H - NH1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NH1'
      bnlocal(isame,3) = 'CT1'
C     --- H - NH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NH1'
      bnlocal(isame,3) = 'CT2'
C     --- H - NH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NH1'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (23.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H - NH2 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NH2'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.50d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.150d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H - NR1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NR1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 127.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.140d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H - NR1 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NR1'
      bnlocal(isame,3) = 'CPH2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.130d0
      bclocal(3) = (15.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H - NR3 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NR3'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.090d0
      bclocal(3) = (15.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H - NR3 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NR3'
      bnlocal(isame,3) = 'CPH2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.00d0
      bclocal(1) = (28.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H - NY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NY'
      bnlocal(isame,3) = 'CA'
C     --- H - NY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'NY'
      bnlocal(isame,3) = 'CPT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H - OH1 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'OH1'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - C - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - C - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3A - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3A'
      bnlocal(isame,3) = 'CN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (29.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1525d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CA - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (41.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CA - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CPT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.00d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.173d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CA - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CY'
C     --- HP - CA - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HP'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CY'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.44d0
      bclocal(1) = (12.70d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CPM - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CPM'
      bnlocal(isame,3) = 'CPA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.40d0
      bclocal(1) = (35.50d0)*(CON_KCALTOK)
      bclocal(2) = 1.7750d0
      bclocal(3) = (14.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CS - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CS'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.50d0
      bclocal(1) = (49.30d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CT2 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CA'
C     --- HA - CT3 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CT2 - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CPB'
C     --- HA - CT3 - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CPB'
C     --- HB - CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'C'
C     --- HB - CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CC'
C     --- HB - CT1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CD'
C     --- HB - CT2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'C'
C     --- HB - CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CC'
C     --- HB - CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CD'
C     --- HB - CT3 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (33.43d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CT2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CPH1'
C     --- HA - CT2 - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CY'
C     --- HA - CT3 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.40d0
      bclocal(1) = (20.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.186d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CY'
      bnlocal(isame,3) = 'CA'
C     --- HP - CY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HP'
      bnlocal(isame,2) = 'CY'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.40d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.255d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CY'
      bnlocal(isame,3) = 'CPT'
C     --- HP - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HP'
      bnlocal(isame,2) = 'CY'
      bnlocal(isame,3) = 'CPT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HB - CP1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'C'
C     --- HB - CP1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CC'
C     --- HB - CP1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HB - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CP2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HB - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT1'
C     --- HB - CT1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT2'
C     --- HB - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.00d0
      bclocal(1) = (36.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HB - CT2 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HB'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (49.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - NC2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NC2'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (40.40d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - NC2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NC2'
      bnlocal(isame,3) = 'CT2'
C     --- HC - NC2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NC2'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - NC2 - HC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NC2'
      bnlocal(isame,3) = 'HC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - NH2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NH2'
      bnlocal(isame,3) = 'CT2'
C     --- HC - NH2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NH2'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.50d0
      bclocal(1) = (39.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - NH2 - HC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NH2'
      bnlocal(isame,3) = 'HC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.074d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - NH3 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NH3'
      bnlocal(isame,3) = 'CT1'
C     --- HC - NH3 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NH3'
      bnlocal(isame,3) = 'CT2'
C     --- HC - NH3 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NH3'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - NH3 - HC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NH3'
      bnlocal(isame,3) = 'HC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.50d0
      bclocal(1) = (51.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - NP - HC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'NP'
      bnlocal(isame,3) = 'HC'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1525d0
      bclocal(3) = (22.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HP - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HP'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.1460d0
      bclocal(3) = (22.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HP - CA - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HP'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CPT'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 130.00d0
      bclocal(1) = (22.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2150d0
      bclocal(3) = (15.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HR1 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HR1'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 130.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.2000d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HR3 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HR3'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 95.00d0
      bclocal(1) = (38.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HS - S - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HS'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 95.00d0
      bclocal(1) = (43.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HS - S - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HS'
      bnlocal(isame,2) = 'S'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.50d0
      bclocal(1) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- N - C - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CP1'
C     --- N - C - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CT1'
C     --- N - C - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CT2'
C     --- N - C - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.20d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- N - CP1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'C'
C     --- N - CP1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CC'
C     --- N - CP1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.80d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- N - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CP2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.00d0
      bclocal(1) = (48.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- N - CP1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'HB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- N - CP3 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'CP3'
      bnlocal(isame,3) = 'CP2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (48.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- N - CP3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'N'
      bnlocal(isame,2) = 'CP3'
      bnlocal(isame,3) = 'HA'
C     --- NH1 - CT1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'HB'
C     --- NH1 - CT2 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (52.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.3642d0
      bclocal(3) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NC2 - C - NC2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NC2'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'NC2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.50d0
      bclocal(1) = (67.70d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NC2 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NC2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.50d0
      bclocal(1) = (51.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NC2 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NC2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA'
C     --- NC2 - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NC2'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
C     --- NH3 - CT1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'HB'
C     --- NH3 - CT2 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HB'
C     --- NP - CP1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NP'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'HB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.50d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH1 - C - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CP1'
C     --- NH1 - C - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CT1'
C     --- NH1 - C - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CT2'
C     --- NH1 - C - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH1 - CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'C'
C     --- NH1 - CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CC'
C     --- NH1 - CT1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CD'
C     --- NH1 - CT2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'C'
C     --- NH1 - CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CC'
C     --- NH1 - CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH1 - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT1'
C     --- NH1 - CT1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT2'
C     --- NH1 - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CT3'
C     --- NH1 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3T - CN1 - NN2U [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3T'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'NN2U'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (51.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH1 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA'
C     --- NH1 - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.50d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH2 - CC - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CP1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 116.50d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.450d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH2 - CC - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT1'
C     --- NH2 - CC - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT2'
C     --- NH2 - CC - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.00d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      bclocal(2) = 1.980d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH2 - CC - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 109.50d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.140d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH2 - CT2 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HB'
C     --- NH2 - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH2'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.00d0
      bclocal(1) = (43.70d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH3 - CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'C'
C     --- NH3 - CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT1'
      bnlocal(isame,3) = 'CC'
C     --- NH3 - CT2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'C'
C     --- NH3 - CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CC'
C     --- NH3 - CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.50d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.101d0
      bclocal(3) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH3 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA'
C     --- NH3 - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NH3'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NP - CP1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NP'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'C'
C     --- NP - CP1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NP'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CC'
C     --- NP - CP1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NP'
      bnlocal(isame,2) = 'CP1'
      bnlocal(isame,3) = 'CD'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.15d0
      bclocal(1) = (51.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NP - CP3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NP'
      bnlocal(isame,2) = 'CP3'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.54d0
      bclocal(1) = (122.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NPH - CPA - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NPH'
      bnlocal(isame,2) = 'CPA'
      bnlocal(isame,3) = 'CPB'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.39d0
      bclocal(1) = (88.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NPH - CPA - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NPH'
      bnlocal(isame,2) = 'CPA'
      bnlocal(isame,3) = 'CPM'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.00d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR1 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR1'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 124.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.140d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR1 - CPH1 - HR3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR1'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'HR3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.50d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.140d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR1 - CPH2 - HR1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR1'
      bnlocal(isame,2) = 'CPH2'
      bnlocal(isame,3) = 'HR1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.00d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR2 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR2'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR2 - CPH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR2'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 120.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.140d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR2 - CPH1 - HR3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR2'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'HR3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.00d0
      bclocal(1) = (25.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.120d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR2 - CPH2 - HR1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR2'
      bnlocal(isame,2) = 'CPH2'
      bnlocal(isame,3) = 'HR1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.50d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR2 - CPH2 - NR1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR2'
      bnlocal(isame,2) = 'CPH2'
      bnlocal(isame,3) = 'NR1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (45.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR3 - CPH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR3'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.00d0
      bclocal(1) = (22.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.180d0
      bclocal(3) = (15.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR3 - CPH1 - HR1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR3'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'HR1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 126.00d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.140d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR3 - CPH1 - HR1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR3'
      bnlocal(isame,2) = 'CPH2'
      bnlocal(isame,3) = 'HR2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.00d0
      bclocal(1) = (120.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.240d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NY - CA - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NY'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CY'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 125.00d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.177d0
      bclocal(3) = (25.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NY - CA - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NY'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'HA'
C     --- NY - CA - HP [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NY'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'HP'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- O - C - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CP1'
C     --- O - CC - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CP1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.00d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- O - C - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CT1'
C     --- O - C - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CT2'
C     --- O - C - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.00d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3A - HN3B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3A'
      bnlocal(isame,3) = 'HN3B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.70d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- O - C - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.50d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- O - C - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'N'
C     --- O - C - NH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'NH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 121.00d0
      bclocal(1) = (15.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.440d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- O - CC - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT1'
C     --- O - CC - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT2'
C     --- O - CC - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- O - CC - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 122.50d0
      bclocal(1) = (75.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.370d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- O - CC - NH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'NH2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OC - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OC'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OC - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OC'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.30d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OC - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OC'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA'
C     --- OC - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OC'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (45.20d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OH1 - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OH1'
      bnlocal(isame,2) = 'CA'
      bnlocal(isame,3) = 'CA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.50d0
      bclocal(1) = (58.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- S - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'S'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT1'
C     --- SM - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'SM'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.50d0
      bclocal(1) = (58.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- S - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'S'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
C     --- S - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'S'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.30d0
      bclocal(1) = (46.10d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- S - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'S'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA'
C     --- S - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'S'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.00d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- SM - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'SM'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA'
C     --- SM - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'SM'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.30d0
      bclocal(1) = (72.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- SM - SM - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'SM'
      bnlocal(isame,2) = 'SM'
      bnlocal(isame,3) = 'CT2'
C     --- SM - SM - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'SM'
      bnlocal(isame,2) = 'SM'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- SS - CS -  CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'SS'
      bnlocal(isame,2) = 'CS'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.30d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- SS - CS - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'SS'
      bnlocal(isame,2) = 'CS'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     *** Nucleic acid angle types (if not listed above) ***
      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.00d0
      bclocal(1) = (42.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2 - CN3C - HN6 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2'
      bnlocal(isame,2) = 'CN3C'
      bnlocal(isame,3) = 'HN6'
C     --- HN6 - CN3C - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN6'
      bnlocal(isame,2) = 'CN3C'
      bnlocal(isame,3) = 'CN3'
C     --- CN3C - CN3 - HN6 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3C'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN6'
C     --- CN3 - CN3C - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3C'
      bnlocal(isame,3) = 'HN3'
C     --- HN3 - CN3C - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN3'
      bnlocal(isame,2) = 'CN3C'
      bnlocal(isame,3) = 'NN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN6 - CN3 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN6'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN8'
C     --- CN8 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.80d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - NN3A - CN4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'NN3A'
      bnlocal(isame,3) = 'CN4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 133.00d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN3A - CN4 - NN3A [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN3A'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'NN3A'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.10d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN3A - CN5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN3A'
      bnlocal(isame,3) = 'CN5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.40d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - CN5 - NN3A [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN3A'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.00d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - CN5 - CN5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'CN5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.70d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - CN2 - NN3A [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN3A'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.70d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - CN5 - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.00d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - CN5 - NN4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN4'
C     --- CN7 - CN7 - ON6B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'ON6B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.60d0
      bclocal(1) = (120.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN4 - CN5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN4'
      bnlocal(isame,3) = 'CN5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.40d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2 - CN4 - NN4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'NN4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.30d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN2 - CN5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.90d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2 - CN5 - NN3A [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN3A'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.00d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - CN5 - NN4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.50d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN3 - CN4 - NN3A [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN3'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'NN3A'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.70d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN3A - CN2 - NN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN3A'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.60d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - CN2 - NN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.50d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - NN1 - HN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'NN1'
      bnlocal(isame,3) = 'HN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (31.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN1 - NN1 - HN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN1'
      bnlocal(isame,2) = 'NN1'
      bnlocal(isame,3) = 'HN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.80d0
      bclocal(1) = (39.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN4 - CN4 - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN4'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'HN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.80d0
      bclocal(1) = (39.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2 - CN4 - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'HN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.40d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - NN2 - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN2 - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 131.10d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - NN2G - CN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2G'
      bnlocal(isame,3) = 'CN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.20d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2G - CN2 - NN3G [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2G'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN3G'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.40d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - NN3G - CN5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'NN3G'
      bnlocal(isame,3) = 'CN5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.90d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5G - CN5 - NN3G [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5G'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN3G'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.60d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - CN5G - CN5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'CN5G'
      bnlocal(isame,3) = 'CN5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.80d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5G - CN1 - NN2G [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5G'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'NN2G'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.60d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5G - CN5 - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5G'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.40d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - CN5G - NN4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'CN5G'
      bnlocal(isame,3) = 'NN4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.80d0
      bclocal(1) = (120.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN4 - CN5G [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN4'
      bnlocal(isame,3) = 'CN5G'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.00d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2B - CN4 - NN4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2B'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'NN4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.20d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN2B - CN5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.50d0
      bclocal(1) = (140.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2B - CN5 - NN3G [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2B'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN3G'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.00d0
      bclocal(1) = (125.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - CN5G - NN4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'CN5G'
      bnlocal(isame,3) = 'NN4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.30d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - NN2G - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2G'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.60d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - NN2G - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'NN2G'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.40d0
      bclocal(1) = (95.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN1 - CN2 - NN2G [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN1'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN2G'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.40d0
      bclocal(1) = (95.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN1 - CN2 - NN3G [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN1'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN3G'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.50d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2G - CN1 - ON1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2G'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'ON1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.70d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5G - CN1 - ON1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5G'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'ON1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.20d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN3 - CN4 - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN3'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.60d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN2B - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.30d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - NN2B - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.10d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - NN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.80d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2 - CN1 - NN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'NN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.10d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - NN3 -CN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN3'
      bnlocal(isame,3) = 'CN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.30d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN2 - NN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.80d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - CN3 - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.90d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3 - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'NN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.20d0
      bclocal(1) = (37.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - NN2 - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.70d0
      bclocal(1) = (37.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - NN2 - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.40d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2 - CN1 - ON1C [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'ON1C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.80d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN3 - CN1 - ON1C [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN3'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'ON1C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.30d0
      bclocal(1) = (81.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN3 - CN2 - NN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN3'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.40d0
      bclocal(1) = (81.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN2 - NN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.10d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.10d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.00d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN3 - CN3 - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN3'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'NN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1T - NN2B - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1T'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN3'
C     --- CN3 - NN2B - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2B - CN1T - NN2U [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2B'
      bnlocal(isame,2) = 'CN1T'
      bnlocal(isame,3) = 'NN2U'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.20d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1T - NN2U - CN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1T'
      bnlocal(isame,2) = 'NN2U'
      bnlocal(isame,3) = 'CN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.60d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2U - CN1 - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2U'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'CN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.60d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - CN3 - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.60d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3 - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (40.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1T - NN2B - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1T'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.00d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - NN2B - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.60d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2B - CN1T - ON1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2B'
      bnlocal(isame,2) = 'CN1T'
      bnlocal(isame,3) = 'ON1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.40d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2U - CN1T - ON1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2U'
      bnlocal(isame,2) = 'CN1T'
      bnlocal(isame,3) = 'ON1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.40d0
      bclocal(1) = (40.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1T - NN2U - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1T'
      bnlocal(isame,2) = 'NN2U'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.40d0
      bclocal(1) = (40.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - NN2U - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2U'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.90d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2U - CN1 - ON1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2U'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'ON1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.50d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN1 - ON1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'ON1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.30d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.30d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN3 - CN3 - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN3'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.70d0
      bclocal(1) = (120.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - CN3T - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'CN3T'
      bnlocal(isame,3) = 'CN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.60d0
      bclocal(1) = (120.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3T - CN3 - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3T'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.60d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3T - CN1 - ON1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3T'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'ON1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.70d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - CN3T - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'CN3T'
      bnlocal(isame,3) = 'CN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.60d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3T - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3T'
      bnlocal(isame,3) = 'CN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.10d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3T - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3T'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.75d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN3 - CN4 - NN3I [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN3'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'NN3I'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.75d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN3 - CN4 - NN2G [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN3'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'NN2G'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.50d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN3I - CN4 - NN2G [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN3I'
      bnlocal(isame,2) = 'CN4'
      bnlocal(isame,3) = 'NN2G'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.00d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5G - CN5 - NN3I [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5G'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN3I'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.40d0
      bclocal(1) = (140.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2B - CN5 - NN3I [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2B'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN3I'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.70d0
      bclocal(1) = (125.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - NN2G - CN4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2G'
      bnlocal(isame,3) = 'CN4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (40.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN2G - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2G'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.40d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN3I - CN5 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN3I'
      bnlocal(isame,3) = 'CN5'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.00d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1T - NN2B - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1T'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.40d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - NN2 - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - NN2 - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.90d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - NN2 - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN9'
C     --- CN5 - NN2B - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN9'
C     --- CN5 - NN2B - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN8'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.80d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN2 - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.90d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN2B - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN9'
C     --- CN4 - NN2B - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN8'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.70d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2B - CN8 - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2B'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'CN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.40d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1T - NN2B - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1T'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.60d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - NN2B - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1 - NN2 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.90d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - NN2 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.10d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - NN2 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.60d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN2 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.50d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN5 - NN2B - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.30d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN4 - NN2B - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN4'
      bnlocal(isame,2) = 'NN2B'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.00d0
      bclocal(1) = (110.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON6B - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2'
C     --- ON6B - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.70d0
      bclocal(1) = (110.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN8 - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2'
C     --- CN8 - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.00d0
      bclocal(1) = (110.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7B - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2'
C     --- CN7C - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7C'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2'
C     --- CN7C - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7C'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2B'
C     --- CN7B - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.00d0
      bclocal(1) = (43.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN7 - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2'
C     --- HN7 - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 113.00d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.460d0
      bclocal(3) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN1A - NN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN1A'
      bnlocal(isame,3) = 'NN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 118.50d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.430d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN1A - ON1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN1A'
      bnlocal(isame,3) = 'ON1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 128.50d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.170d0
      bclocal(3) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN1 - CN1A - ON1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NN1'
      bnlocal(isame,2) = 'CN1A'
      bnlocal(isame,3) = 'ON1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.40d0
      bclocal(1) = (110.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3D - CN2 - NN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3D'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.30d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3D - CN2 - NN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3D'
      bnlocal(isame,2) = 'CN2'
      bnlocal(isame,3) = 'NN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.80d0
      bclocal(1) = (5.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1A - CN3 - CN3C [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1A'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN3C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.20d0
      bclocal(1) = (125.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1A - CN3 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1A'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN8'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.80d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - CN3D - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'CN3D'
      bnlocal(isame,3) = 'CN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN2 - CN3D - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN2'
      bnlocal(isame,2) = 'CN3D'
      bnlocal(isame,3) = 'CN9'
C     --- CN3A - CN3 - CN3B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3A'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN3B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.20d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3A - CN3 - CN1A [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3A'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN1A'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 131.80d0
      bclocal(1) = (10.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3B - CN3 - CN1A [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3B'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN1A'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.50d0
      bclocal(1) = (53.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN8'
C     --- CN3B - CN3 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3B'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN8'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.00d0
      bclocal(1) = (43.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3C - CN3 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3C'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'CN8'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.20d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3D - CN9 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3D'
      bnlocal(isame,3) = 'CN9'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3D - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3D'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN3'
C     --- CN3B - CN3 - HN3B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3B'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN3B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3A - CN3 - HN3B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3A'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN3B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.00d0
      bclocal(1) = (42.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3C - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3C'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'HN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3B - HN3B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3B'
      bnlocal(isame,3) = 'HN3B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (120.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN3B - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN3B'
      bnlocal(isame,3) = 'NN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.20d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3D - CN3 - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3D'
      bnlocal(isame,2) = 'CN3'
      bnlocal(isame,3) = 'NN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.50d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN3B - CN3B - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN3B'
      bnlocal(isame,2) = 'CN3B'
      bnlocal(isame,3) = 'NN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (125.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN8 - CN3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'CN3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.10d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3 - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'HN8'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.10d0
      bclocal(1) = (33.43d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN8 - CN8 - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'NN2'
C     --- HN9 - CN9 - NN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN9'
      bnlocal(isame,2) = 'CN9'
      bnlocal(isame,3) = 'NN2'
C     --- HN9 - CN9 - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN9'
      bnlocal(isame,2) = 'CN9'
      bnlocal(isame,3) = 'NN2B'
C     --- HN8 - CN8 - NN2B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'NN2B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN1A - NN1 - HN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1A'
      bnlocal(isame,2) = 'NN1'
      bnlocal(isame,3) = 'HN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3B - NN2 - CN3B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3B'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN3B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.00d0
      bclocal(1) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3C - NN2 - CN3C [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3C'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN3C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.70d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3B - NN2 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3B'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN7B'
C     --- CN3C - NN2 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3C'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.40d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3B - NN2 - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3B'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.00d0
      bclocal(1) = (39.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN3C - NN2 - HN2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN3C'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'HN2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.00d0
      bclocal(1) = (39.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN1 - NN5 - HN1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN1'
      bnlocal(isame,2) = 'NN5'
      bnlocal(isame,3) = 'HN1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 140.00d0
      bclocal(1) = (15.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.800d0
      bclocal(3) = (-40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- P - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'P'
      bnlocal(isame,2) = 'ON2'
      bnlocal(isame,3) = 'P'
C     --- P2 - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'P2'
      bnlocal(isame,2) = 'ON2'
      bnlocal(isame,3) = 'P'
C     --- P2 - ON2 - P2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'P2'
      bnlocal(isame,2) = 'ON2'
      bnlocal(isame,3) = 'P2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 115.00d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.350d0
      bclocal(3) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN4 - ON4 - P [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN4'
      bnlocal(isame,2) = 'ON4'
      bnlocal(isame,3) = 'P'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.00d0
      bclocal(1) = (98.90d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON4 - P - ON4 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON4'
      bnlocal(isame,2) = 'P'
      bnlocal(isame,3) = 'ON4'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.60d0
      bclocal(1) = (88.90d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON2 - P2 - ON3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON2'
      bnlocal(isame,2) = 'P2'
      bnlocal(isame,3) = 'ON3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (104.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON3 - P2 - ON3 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON3'
      bnlocal(isame,2) = 'P2'
      bnlocal(isame,3) = 'ON3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 102.00d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON6 - CN7B - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CN8'
C     --- ON6 - CN7B - CN7C [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CN7C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 100.00d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7B - CN8 - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'CN7'
C     --- CN7B - CN7C - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7C'
      bnlocal(isame,3) = 'CN7'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 102.00d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.561d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN8 - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7'
C     --- CN7C - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7C'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.00d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7 - CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'ON6'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.24d0
      bclocal(1) = (45.20d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN7 - CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'ON6'
C     --- HN7 - CN7B - ON6 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'ON6'
C     --- HN7 - CN7B - ON6B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'ON6B'
C     --- HN7 - CN7 - ON6B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'ON6B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN7 - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (34.53d0)*(CON_KCALTOK)
      bclocal(2) = 2.179d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN8 - CN8 - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'CN7'
C     --- HN8 - CN8 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8'
      bnlocal(isame,3) = 'CN8'
C     --- HN7 - CN7C - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7C'
      bnlocal(isame,3) = 'CN7'
C     --- HN8 - CN8B - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN8'
      bnlocal(isame,2) = 'CN8B'
      bnlocal(isame,3) = 'CN7'
C     --- HN7 - CN7B - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CN7'
C     --- HN7 - CN7 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.00d0
      bclocal(1) = (110.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7 - ON6B - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'ON6B'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.00d0
      bclocal(1) = (90.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON6B - CN7B - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 96.00d0
      bclocal(1) = (110.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7B - CN7B - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CN7'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 100.00d0
      bclocal(1) = (60.0d0)*(CON_KCALTOK)
      bclocal(2) = 2.561d0
      bclocal(3) = (8.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7B - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 110.10d0
      bclocal(1) = (33.4d0)*(CON_KCALTOK)
      bclocal(2) = 2.179d0
      bclocal(3) = (22.53d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7B - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'HN7'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.20d0
      bclocal(1) = (90.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON6 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN8B'
C     --- ON6B - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6B'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN8B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7 - CN8B - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'CN8B'
      bnlocal(isame,3) = 'ON2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.70d0
      bclocal(1) = (115.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7 - CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'ON2'
C     --- CN7B - CN7B - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'ON2'
C     --- CN8 - CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'ON2'
C     --- CN7 - CN7B - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'ON2'
C     --- CN7C - CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7C'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'ON2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.00d0
      bclocal(1) = (75.7d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON5 - CN7 - CN8 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON5'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN8'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.00d0
      bclocal(1) = (57.5d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN5 - ON5 - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN5'
      bnlocal(isame,2) = 'ON5'
      bnlocal(isame,3) = 'CN7'
C     --- HN5 - ON5 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN5'
      bnlocal(isame,2) = 'ON5'
      bnlocal(isame,3) = 'CN7B'
C     --- HN5 - ON5 - CN7C [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN5'
      bnlocal(isame,2) = 'ON5'
      bnlocal(isame,3) = 'CN7C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.00d0
      bclocal(1) = (90.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON2 - CN7 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON2'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7B'
C     --- ON5 - CN7 - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON5'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.40d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON5 - CN7B - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON5'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CN7B'
C     --- ON5 - CN7C - CN7B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON5'
      bnlocal(isame,2) = 'CN7C'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (90.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON5 - CN7B - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON5'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CN7'
C     --- ON5 - CN7C - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON5'
      bnlocal(isame,2) = 'CN7C'
      bnlocal(isame,3) = 'CN7'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (140.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON6 - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NR1'
C     --- ON6B - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NR1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.70d0
      bclocal(1) = (140.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN8 - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN8'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NR1'
C     --- CN7B - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NR1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.00d0
      bclocal(1) = (130.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7B - NR1 - CPH2 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'NR1'
      bnlocal(isame,3) = 'CPH2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.00d0
      bclocal(1) = (130.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7B - NR1 - CPH1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7B'
      bnlocal(isame,2) = 'NR1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.30d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN7 - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NR1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.00d0
      bclocal(1) = (130.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CPH2 - NR1 - CPH1 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CPH2'
      bnlocal(isame,2) = 'NR1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.25d0
      bclocal(1) = (120.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- ON6 - CN7B - CN7 [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'ON6'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CN7'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.00d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN7 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CN7'
      bnlocal(isame,2) = 'CN7'
      bnlocal(isame,3) = 'CN8B'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.05d0
      bclocal(1) = (96.150d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- FE - NPH - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'FE'
      bnlocal(isame,2) = 'NPH'
      bnlocal(isame,3) = 'CPA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 133.0d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- FE - NR2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'FE'
      bnlocal(isame,2) = 'NR2'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.0d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- FE - NR2 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'FE'
      bnlocal(isame,2) = 'NR2'
      bnlocal(isame,3) = 'CPH2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 180.0d0
      bclocal(1) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CPM - FE [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CPM'
      bnlocal(isame,3) = 'FE'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 90.0d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NPH - FE - CM [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NPH'
      bnlocal(isame,2) = 'FE'
      bnlocal(isame,3) = 'CM'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 45.0d0
      bclocal(1) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NPH - FE - CPM [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NPH'
      bnlocal(isame,2) = 'FE'
      bnlocal(isame,3) = 'CPM'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 90.0d0
      bclocal(1) = (14.390d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NPH - FE - NPH [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NPH'
      bnlocal(isame,2) = 'FE'
      bnlocal(isame,3) = 'NPH'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 180.0d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR2 - FE - CM [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR2'
      bnlocal(isame,2) = 'FE'
      bnlocal(isame,3) = 'CM'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 90.0d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR2 - FE - NPH [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'NR2'
      bnlocal(isame,2) = 'FE'
      bnlocal(isame,3) = 'NPH'
c     --- set nbsame
      call twh_set_bendnames(type,nbslocal,bnlocal)
      nbslocal = isame

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 180.0d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OM - CM - FE [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OM'
      bnlocal(isame,2) = 'CM'
      bnlocal(isame,3) = 'FE'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 90.0d0
      bclocal(1) = (5.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OM - FE - NPH [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OM'
      bnlocal(isame,2) = 'FE'
      bnlocal(isame,3) = 'NPH'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 180.0d0
      bclocal(1) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OM - OM - FE [par_all27_prot_lipid.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'OM'
      bnlocal(isame,2) = 'OM'
      bnlocal(isame,3) = 'FE'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 121.70d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 0.0d0
      bclocal(3) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- O - C - H
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'H'
C     --- O - C - HA (ethanolamine cap) References for this group are
c     --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c         `High-resolution polypeptide structure in a lamerllar 
c          phase lipid environment from solid state NMR derived 
c          orientational constraints'
c         Structure, 5, 1655-1669 (1997).
c     --- T.B. Woolf personal communication 08-16-2001 (M.G. Martin)
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 111.00d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      bclocal(2) = 1.9800d0
      bclocal(3) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NH2 - CC - HA
      isame = isame + 1
      bnlocal(isame,1) = 'NH2'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'HA'
C     --- NH1 - C - HA (CHO parameters for formamide cap on gramicidin)
c     --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c         `High-resolution polypeptide structure in a lamellar 
c          phase lipid environment from solid state NMR derived 
c          orientational constraints'
c         Structure, 5, 1655-1669 (1997).
c     --- T.B. Woolf personal communication 08-16-2001 (M.G. Martin)
      isame = isame + 1
      bnlocal(isame,1) = 'NH1'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return
      return
      end

      subroutine twh_ffcharmm27_dihedral(lfailure,ffname,ndihedtypes)
c     ******************************************************************
c     * Charmm27 standard  dihedral                                    *
c     *                                                                *
c     * split from original ffcharmm27 03-05-2006 by M.G. Martin       *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer ndihedtypes
c     --- local variables
      integer type,isame
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision onepi

      call twh_constant_pi(GLB_GET,onepi)

c     --- initialize type count
      type = 0

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
      isame = 0
c     --- X - CTL1 - CTL1 - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL1'
      tnlocal(isame,3) = 'CTL1'
      tnlocal(isame,4) = 'wild'
c     --- X - CTL1 - CTL2 - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'wild'
c     --- X - CTL1 - CTL3 - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL1'
      tnlocal(isame,3) = 'CTL3'
      tnlocal(isame,4) = 'wild'
c     --- C - CT1 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'C'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- C - CT2 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'C'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- CC - CT1 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CC'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- CC - CT2 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CC'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- CD - CT1 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CD'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- CD - CT2 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CD'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)
      
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
      isame = 0

c     --- OBL - CL - CTL2 - HAL2 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OBL'
      tnlocal(isame,2) = 'CL'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'HAL2'
c     --- OBL - CL - CTL3 - HAL3 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OBL'
      tnlocal(isame,2) = 'CL'
      tnlocal(isame,3) = 'CTL3'
      tnlocal(isame,4) = 'HAL3'
c     --- OSL - CL - CTL2 - HAL2 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OSL'
      tnlocal(isame,2) = 'CL'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'HAL2'
c     --- OSL - CL - CTL3 - HAL3 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OSL'
      tnlocal(isame,2) = 'CL'
      tnlocal(isame,3) = 'CTL3'
      tnlocal(isame,4) = 'HAL3'
c     --- HAL2 - CTL2 - CL - OHL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HAL2'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CL'
      tnlocal(isame,4) = 'OHL'
c     --- HAL3 - CTL3 - CL - OHL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HAL3'
      tnlocal(isame,2) = 'CTL3'
      tnlocal(isame,3) = 'CL'
      tnlocal(isame,4) = 'OHL'
c     --- X - CTL1 - OSL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL1'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'wild'
c     --- X - CTL2 - OSL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'wild'
c     --- X - CTL3 - OSL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL3'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'wild'
c     --- HAL3 - CTL3 - OSL - SL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HAL3'
      tnlocal(isame,2) = 'CTL3'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'SL'
c     --- CTL2 - OSL - SL - O2L [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL2'
      tnlocal(isame,2) = 'OSL'
      tnlocal(isame,3) = 'SL'
      tnlocal(isame,4) = 'O2L'
c     --- CTL3 - OSL - SL - O2L [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL3'
      tnlocal(isame,2) = 'OSL'
      tnlocal(isame,3) = 'SL'
      tnlocal(isame,4) = 'O2L'
c     --- HEL1 - CEL1 - CTL2 - HAL2 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'HAL2'
c     --- HEL1 - CEL1 - CTL3 - HAL3 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL3'
      tnlocal(isame,4) = 'HAL3'
c     --- CEL1 - CTL2 - CEL1 - HEL1 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL1'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CEL1'
      tnlocal(isame,4) = 'HEL1'
c     --- CD - CP1 - N - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CD'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'C'
c     --- CP2 - CP3 - N - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP2'
      tnlocal(isame,2) = 'CP3'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'C'
c     --- H - NH1 - CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'C'
c     --- H - NH1 - CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CC'
c     --- H - NH1 - CT1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CD'
c     --- H - NH1 - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT1'
c     --- H - NH1 - CT1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT2'
c     --- H - NH1 - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT3'
c     --- H - NH1 - CT2 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'C'
c     --- H - NH1 - CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CC'
c     --- H - NH1 - CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CD'
c     --- H - NH1 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- H - NH1 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT3'
c     --- HA - CP3 - N - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CP3'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'C'
c     --- HA - CT2 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- HA - CT2 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- HA - CT2 - NH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'H'
c     --- HA - CT3 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- HA - CT3 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- HA - CT3 - NH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'H'
c     --- HB - CT1 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HB'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- HB - CT1 - NH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HB'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'H'
c     --- HB - CT2 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HB'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- HB - CT2 - NH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HB'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'H'
c     --- HB - CT3 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HB'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- HB - CT3 - NH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HB'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'H'
c     --- HR2 - CPH2 - NR3 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR2'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR3'
      tnlocal(isame,4) = 'H'
c     --- N - C - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT1'
c     --- N - C - CT1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT2'
c     --- N - C - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT3'
c     --- N - C - CT1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'HB'
c     --- N - C - CT2 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HB'
c     --- N - C - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- NH1 - C - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT1'
c     --- NH1 - C - CT1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT2'
c     --- NH1 - C - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT3'
c     --- NH1 - C - CT1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'HB'
c     --- NH1 - C - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- NH1 - C - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- NH1 - C - CT2 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HB'
c     --- NH1 - C - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- NH2 - CC - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH2'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- O - C - CT1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'HB'
c     --- O - C - CT1 - NH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'NH1'
c     --- O - C - CT1 - NH3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'NH3'
c     --- O - C - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- O - C - CT2 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HB'
c     --- O - C - CT2 - NH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'NH1'
c     --- O - C - CT2 - NH3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'NH3'
c     --- O - C - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- O - CC - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- X - CP1 - C - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'wild'
c     --- X - CP1 - CC - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'wild'
c     --- X - CP1 - CD - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'CD'
      tnlocal(isame,4) = 'wild'
c     --- X - CPA - CPB - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CPA'
      tnlocal(isame,3) = 'CPB'
      tnlocal(isame,4) = 'wild'
c     --- X - CPA - CPM - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CPA'
      tnlocal(isame,3) = 'CPM'
      tnlocal(isame,4) = 'wild'
c     --- X - CPB - CPB - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CPB'
      tnlocal(isame,3) = 'CPB'
      tnlocal(isame,4) = 'wild'
c     --- X - CPB - CT2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CPB'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'wild'
c     --- X - CPB - CT3 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CPB'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'wild'
c     --- X - CPT - CPT - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'wild'
c     --- X - CT1 - CD - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'CD'
      tnlocal(isame,4) = 'wild'
c     --- X - CT2 - CA - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'wild'
c     --- X - CT2 - CD - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CD'
      tnlocal(isame,4) = 'wild'
c     --- X - CT2 - NC2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NC2'
      tnlocal(isame,4) = 'wild'
c     --- X - CT3 - CA - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'wild'
c     --- X - CT3 - CD - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'CD'
      tnlocal(isame,4) = 'wild'
c     --- X - CT3 - NC2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'NC2'
      tnlocal(isame,4) = 'wild'
c     --- X - NPH - CPA - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'NPH'
      tnlocal(isame,3) = 'CPA'
      tnlocal(isame,4) = 'wild'
c     --- P - ON2 - CN7 - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'HN7'
c     --- P - ON2 - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- P - ON2 - CN8B - HN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'HN8'
c     --- P - ON2 - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'HN8'
c     --- P - ON2 - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- ON1 - CN1T - NN2B - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1'
      tnlocal(isame,2) = 'CN1T'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'HN2'
c     --- ON1 - CN1 - NN2U - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2U'
      tnlocal(isame,4) = 'HN2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     --- null2 : null torsions split into two groups
c     --- in order to keep maxsamet at 90
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
      isame = 0
c     --- ON1 - CN1T - NN2U - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1'
      tnlocal(isame,2) = 'CN1T'
      tnlocal(isame,3) = 'NN2U'
      tnlocal(isame,4) = 'HN2'
c     --- NN1 - CN2 - CN3D - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN1'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN3D'
      tnlocal(isame,4) = 'CN3'
c     --- NN1 - CN2 - CN3D - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN1'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN3D'
      tnlocal(isame,4) = 'CN9'
c     --- NN3 - CN2 - CN3D - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN3D'
      tnlocal(isame,4) = 'CN9'
c     --- CN3 - CN3D - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN3D'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- NN4 - CN5 - CN2 - NN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN4'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN2'
      tnlocal(isame,4) = 'NN1'
c     --- HN2 - NN2 - CN4 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN2'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'HN3'
c     --- ON1 - CN1 - CN5G - NN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN5G'
      tnlocal(isame,4) = 'NN4'
c     --- ON1 - CN1 - NN2G - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'HN2'
c     --- NN1 - CN2 - NN2G - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN1'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'HN2'
c     --- HN3 - CN4 - NN2B - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN3'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'HN2'
c     --- X - CN5 - CN5 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'wild'
c     --- X - CN5G - CN5 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN5G'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'wild'
c     --- CN3 - NN2 - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN4 - NN2 - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN3 - NN2B - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN4 - NN2B - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN4 - NN2B - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'HN8'
c     --- CN4 - NN2B - CN8 - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN9'
c     --- CN7 - CN7B - NN2 - CN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3B'
c     --- CN7B - CN7B - NN2 - CN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3B'
c     --- ON6 - CN7B - NN2 - CN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3B'
c     --- ON6B - CN7B - NN2 - CN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3B'
c     --- CN7 - CN7B - NN2 - CN3C [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3C'
c     --- CN7B - CN7B - NN2 - CN3C [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3C'
c     --- ON6 - CN7B - NN2 - CN3C [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3C'
c     --- ON6B - CN7B - NN2 - CN3C [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3C'
c     --- CN7 - CN7 - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- CN7 - CN7B - NN2 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN4'
c     --- ON2 - CN7 - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- ON2 - CN7B - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- ON5 - CN7B - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- ON5 - CN7 - CN7 - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON5'
c     --- ON5 - CN7 - CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON2'
c     --- CN7B - NR1 - CPH1 - HR3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'HR3'
c     --- CN7B - NR1 - CPH1 - CPH1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- CN7B - NR1 - CPH2 - HR1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH2'
      tnlocal(isame,4) = 'HR1'
c     --- CN7B - NR1 - CPH2 - HR2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH2'
      tnlocal(isame,4) = 'HR2'
c     --- HN7 - CN7B - NR1 - CPH1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH1'
c     --- CN8 - CN8 - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NR1'
c     --- CN7 - CN7B - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NR1'
c     --- CN8 - ON6 - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'ON6'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NR1'
c     --- HN8 - CN8 - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NR1'
c     --- HN7 - CN7B - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NR1'
c     --- NR1 - CN7B - CN7B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON5'
c     --- CN7B - ON6 - CN7 - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'ON6'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN9'
c     --- CN7 - ON6 - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'ON6'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NR1'
c     --- CN7 - ON6B - CN7B - NR1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'ON6B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NR1'
c     --- NR1 - CN7B - CN8 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN7'
c     --- CN8 - CN7B - NN2 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN4'
c     --- CN7B - CN7B - NN2 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN4'
c     --- CN7C - CN7B - NN2 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN4'
c     --- HN7 - CN7B - NN2 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN5'
c     --- ON6 - CN7B - NN2 - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN1'
c     --- ON6B - CN7B - NN2 - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN1'
c     --- CN8 - CN7B - NN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3'
c     --- CN7B - CN7B - NN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3'
c     --- CN7C - CN7B - NN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3'
c     --- HN7 - CN7B - NN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN1'
c     --- ON6 - CN7B - NN2B - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN1'
c     --- ON6B - CN7B - NN2B - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN1'
c     --- CN8 - CN7B - NN2B - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN3'
c     --- CN7B - CN7B - NN2B - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN3'
c     --- CN7C - CN7B - NN2B - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN3'
c     --- HN7 - CN7B - NN2B - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN3'
c     --- CN8 - CN7B - NN2B - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN5'
c     --- CN8 - CN7B - NN2B - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN4'
c     --- CN7B - CN7B - NN2B - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN5'
c     --- CN7B - CN7B - NN2B - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN4'
c     --- CN7C - CN7B - NN2B - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN5'
c     --- CN7C - CN7B - NN2B - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN4'
c     --- HN7 - CN7B - NN2B - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN5'
c     --- CN7 - ON6 - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'ON6'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- CN7 - ON6 - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'ON6'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2B'
c     --- CN7 - ON6B - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'ON6B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- CN7 - ON6B - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'ON6B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2B'
c     --- CN7 - CN8 - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2B'
c     --- CN7 - CN8 - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- HN8 - CN8 - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- HN8 - CN8 - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2B'
c     --- CN7 - CN7B - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2B'
c     --- CN7 - CN7B - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- CN7 - CN7C - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2B'
c     --- CN7 - CN7C - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- HN7 - CN7B - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2B'
c     --- HN7 - CN7B - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- HN7 - CN7C - CN7B - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2B'
c     --- HN7 - CN7C - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- HN8 - CN8B - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8B'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- HN7 - CN7 - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     --- null3 : null torsions split into third group
c     --- in order to keep maxsamet at 90
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
      isame = 0
c     --- CN8 - CN7 - CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6'
c     --- CN7B - CN7B - ON6B - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'ON6B'
      tnlocal(isame,4) = 'CN7'
c     --- CN7B - CN7 - CN7 - ON6B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6B'
c     --- CN7C - CN7 - CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6'
c     --- CN7B - CN7B - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- CN7 - CN7 - ON6B - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON6B'
      tnlocal(isame,4) = 'CN7B'
c     --- HN7 - CN7 - CN7 - CN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8'
c     --- HN7 - CN7 - CN7 - CN7C [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7C'
c     --- NN2 - CN7B - CN7B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON5'
c     --- NN2B - CN7B - CN7B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON5'
c     --- ON5 - CN7B - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- HN7 - CN7B - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- CN7 - CN7 - CN7B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON5'
c     --- ON6B - CN7B - CN7B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON5'
c     --- ON5 - CN7B - CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON2'
c     --- ON5 - CN7 - CN7B - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON2'
c     --- ON5 - CN7B - CN7 - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON5'
c     --- HN7 - CN7B - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- CN7B - CN7C - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- CN7 - CN7C - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- ON6 - CN7B - CN7C - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'ON5'
c     --- HN7 - CN7B - CN7C - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'ON5'
c     --- CN7 - CN7 - CN7C - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'ON5'
c     --- HN7 - CN7C - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- ON5 - CN7C - CN7 - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'HN7'
c     --- ON5 - CN7C - CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON2'
c     --- CN7 - ON6 - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'ON6'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- CN1T - NN2B - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1T'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- HN7 - CN7B - ON6B - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'ON6B'
      tnlocal(isame,4) = 'CN7'
c     --- wild - FE - NPH - wild [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'FE'
      tnlocal(isame,3) = 'NPH'
      tnlocal(isame,4) = 'wild'
c     --- wild - FE - OM - wild [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'FE'
      tnlocal(isame,3) = 'OM'
      tnlocal(isame,4) = 'wild'
c     --- wild - NPH - CPA - wild [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'NPH'
      tnlocal(isame,3) = 'CPA'
      tnlocal(isame,4) = 'wild'
c     --- HA - C - NH1 - CT1 (formamide cap for gramicidin)
c     --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c         `High-resolution polypeptide structure in a lamerllar 
c          phase lipid environment from solid state NMR derived 
c          orientational constraints'
c         Structure, 5, 1655-1669 (1997).
c     --- T.B. Woolf personal communication 08-16-2001 (M.G. Martin)
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CTL2 - NH3L - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'NH3L'
      tnlocal(isame,4) = 'wild'
c     --- CP2 - CP3 - N - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP2'
      tnlocal(isame,2) = 'CP3'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP1'
c     --- CP3 - N - CP1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP3'
      tnlocal(isame,2) = 'N'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'C'
c     --- CP3 - N - CP1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP3'
      tnlocal(isame,2) = 'N'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CC'
c     --- CP3 - N - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP3'
      tnlocal(isame,2) = 'N'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CP2'
c     --- HA - CP3 - N - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CP3'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP1'
c     --- HB - CP1 - N - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HB'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP3'
c     --- X - CT1 - NH3 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'NH3'
      tnlocal(isame,4) = 'wild'
c     --- X - CT2 - NH3 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NH3'
      tnlocal(isame,4) = 'wild'
c     --- ON3 - P - ON2 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON3'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN7'
c     --- ON3 - P - ON2 - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON3'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN7B'
c     --- ON3 - P - ON2 - CN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON3'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN8'
c     --- ON3 - P - ON2 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON3'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN8B'
c     --- ON3 - P - ON2 - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON3'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN9'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CTL2 - CTL3 - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL3'
      tnlocal(isame,4) = 'wild'
c     --- HA - CT3 - CS - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'CS'
      tnlocal(isame,4) = 'HA'
c     --- OC - CC - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OC'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CP2'
c     --- OC - CC - CP1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OC'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'HB'
c     --- OC - CC - CP1 - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OC'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'N'
c     --- OC - CC - CP1 - NP [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OC'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'NP'
c     --- X - CP2 - CP2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CP2'
      tnlocal(isame,3) = 'CP2'
      tnlocal(isame,4) = 'wild'
c     --- X - CT2 - CT3 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CTL2 - CTL2 - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'wild'
c     --- NR1 - CPH1 - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT1'
c     --- NR1 - CPH1 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- NR1 - CPH1 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT3'
c     --- NR1 - CPH1 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- NR1 - CPH1 - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- NR2 - CPH1 - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT1'
c     --- NR2 - CPH1 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- NR2 - CPH1 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT3'
c     --- NR2 - CPH1 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- NR2 - CPH1 - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- NR3 - CPH1 - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT1'
c     --- NR3 - CPH1 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- NR3 - CPH1 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT3'
c     --- NR3 - CPH1 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- NR3 - CPH1 - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- CN1 - NN2 - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN5 - NN2 - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN1 - NN2B - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN1T - NN2B - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1T'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN5 - NN2B - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN5 - NN2B - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'HN8'
c     --- CN5 - NN2B - CN8 - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN9'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- OBL - CL - OSL - CTL1 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OBL'
      tnlocal(isame,2) = 'CL'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CTL1'
c     --- OBL - CL - OSL - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OBL'
      tnlocal(isame,2) = 'CL'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CTL2'
c     --- OBL - CL - OSL - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OBL'
      tnlocal(isame,2) = 'CL'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CTL3'
c     --- OB - CD - OS - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OB'
      tnlocal(isame,2) = 'CD'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'CT2'
c     --- OB - CD - OS - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OB'
      tnlocal(isame,2) = 'CD'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'CT3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CL - OSL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CL'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'wild'
c     --- X - CL - OHL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CL'
      tnlocal(isame,3) = 'OHL'
      tnlocal(isame,4) = 'wild'
c     --- X - CD - OH1 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CD'
      tnlocal(isame,3) = 'OH1'
      tnlocal(isame,4) = 'wild'
c     --- X - CD - OS - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CD'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CTL1 - OHL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL1'
      tnlocal(isame,3) = 'OHL'
      tnlocal(isame,4) = 'wild'
c     --- X - CTL2 - OHL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'OHL'
      tnlocal(isame,4) = 'wild'
c     --- X - CTL3 - OHL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL3'
      tnlocal(isame,3) = 'OHL'
      tnlocal(isame,4) = 'wild'
c     --- X - CP1 - CP2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'CP2'
      tnlocal(isame,4) = 'wild'
c     --- X - CP3 - CP2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CP3'
      tnlocal(isame,3) = 'CP2'
      tnlocal(isame,4) = 'wild'
c     --- X - CT1 - OH1 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'OH1'
      tnlocal(isame,4) = 'wild'
c     --- X - CT2 - OH1 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'OH1'
      tnlocal(isame,4) = 'wild'
c     --- X - CT3 - OH1 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'OH1'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CTL2 - CL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CL'
      tnlocal(isame,4) = 'wild'
c     --- X - CTL3 - CL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL3'
      tnlocal(isame,3) = 'CL'
      tnlocal(isame,4) = 'wild'
c     --- X - CT1 - CC - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'wild'
c     --- X - CT2 - CC - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'wild'
c     --- X - CT3 - CC - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(1) = (1.20d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.10d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.10d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      tclocal(10) = (0.10d0)*(CON_KCALTOK)
      tclocal(11) = 3.0d0
      tclocal(12) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- OSL - PL - OSL - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OSL'
      tnlocal(isame,2) = 'PL'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CTL2'
c     --- OSL - PL - OSL - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OSL'
      tnlocal(isame,2) = 'PL'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CTL3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.95d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.50d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- OHL - PL - OSL - CTL2 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OHL'
      tnlocal(isame,2) = 'PL'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CTL2'
c     --- OHL - PL - OSL - CTL3 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OHL'
      tnlocal(isame,2) = 'PL'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CTL3'
c     --- ON4 - P - ON2 - CN7 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON4'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN7'
c     --- ON4 - P - ON2 - CN7B [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON4'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN7B'
c     --- ON4 - P - ON2 - CN8 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON4'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN8'
c     --- ON4 - P - ON2 - CN8B [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON4'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN8B'
c     --- ON4 - P - ON2 - CN9 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON4'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN9'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.30d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- X - OHL - PL - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'OHL'
      tnlocal(isame,3) = 'PL'
      tnlocal(isame,4) = 'wild'
c     --- X - ON4 - P - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'ON4'
      tnlocal(isame,3) = 'P'
      tnlocal(isame,4) = 'wild'
c     --- CN8 - CN7B - NN2 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN5'
c     --- CN7B - CN7B - NN2 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN5'
c     --- CN7C - CN7B - NN2 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.70d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL3 - CTL2 - OSL - CL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL3'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CL'
c     --- CTL2 - CTL2 - OSL - CL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL2'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CL'
c     --- CTL3 - CTL1 - OSL - CL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL3'
      tnlocal(isame,2) = 'CTL1'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CL'
c     --- CTL2 - CTL1 - OSL - CL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL2'
      tnlocal(isame,2) = 'CTL1'
      tnlocal(isame,3) = 'OSL'
      tnlocal(isame,4) = 'CL'
c     --- NH3L - CTL2 - CTL2 - OHL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH3L'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'OHL'
c     --- NH3L - CTL2 - CTL2 - OSL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH3L'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'OSL'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (-0.15d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.53d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL1 - CTL2 - CL - OSL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL1'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CL'
      tnlocal(isame,4) = 'OSL'
c     --- CTL2 - CTL2 - CL - OSL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL2'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CL'
      tnlocal(isame,4) = 'OSL'
c     --- CTL3 - CTL2 - CL - OSL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL3'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CL'
      tnlocal(isame,4) = 'OSL'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.26d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- X - CTL2 - NTL - X  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'NTL'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.23d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- X - CTL5 - NTL - X  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL5'
      tnlocal(isame,3) = 'NTL'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (4.3d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (-0.4d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- NTL - CTL2 - CTL2 - OHL  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NTL'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'OHL'
c     --- NTL - CTL2 - CTL2 - OSL  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NTL'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'OSL'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1525d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- X - CTL3 - CTL3 - X  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CTL3'
      tnlocal(isame,3) = 'CTL3'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.10d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.15d0)*(CON_KCALTOK)
      tclocal(5) = 4.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.10d0)*(CON_KCALTOK)
      tclocal(8) = 6.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL3 - CTL2 - CTL2 - CTL3  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL3'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CTL3'
c     --- CTL2 - CTL2 - CTL2 - CTL3  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL2'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CTL3'
c     --- CTL2 - CTL2 - CTL2 - CTL2  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL2'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CTL2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HEL1 - CEL1 - CEL1 - HEL1  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CEL1'
      tnlocal(isame,4) = 'HEL1'
c     --- CTL3 - CEL1 - CEL1 - HEL1  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL3'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CEL1'
      tnlocal(isame,4) = 'HEL1'
c     --- H - NR1 - CPH1 - CPH1  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- H - NR1 - CPH1 - CT2  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT2'
c     --- H - NR1 - CPH1 - CT3  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT3'
c     --- HA - CA - NY - H  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'NY'
      tnlocal(isame,4) = 'H'
c     --- HP - NY - CA - HP  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'NY'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'HP'
c     --- HR1 - CPH1 - CPH1 - CT2  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT2'
c     --- HR1 - CPH1 - CPH1 - CT3  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT3'
c     --- HR1 - CPH1 - CPH1 - HR1  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'HR1'
c     --- HR1 - CPH2 - NR1 - H  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR1'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'H'
c     --- HR3 - CPH1 - NR1 - H  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'H'
c     --- NR2 - CPH2 - NR1 - H  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'H'
c     --- NN3 - CN2 - NN1 - HN1  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN1'
      tnlocal(isame,4) = 'HN1'
c     --- CN3 - CN2 - NN1 - HN1  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN1'
      tnlocal(isame,4) = 'HN1'
c     --- CN3 - CN3D - CN2 - NN3  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN3D'
      tnlocal(isame,3) = 'CN2'
      tnlocal(isame,4) = 'NN3'
c     --- CN3D - CN3 - NN2 - CN1  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3D'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN1'
c     --- X - CN4 - NN3I - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN3I'
      tnlocal(isame,4) = 'wild'
c     --- X - CN1 - NN3 - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN3'
      tnlocal(isame,4) = 'wild'
c     --- X - CN3 - NN2 - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'wild'
c     --- X - CN3 - NN2B - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'wild'
c     --- X - CN3 - CN3 - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'wild'
c     --- X - CN3 - CN3T - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3T'
      tnlocal(isame,4) = 'wild'
c     --- X - CN1 - CN3 - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'wild'
c     --- X - CN1 - CN3T - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN3T'
      tnlocal(isame,4) = 'wild'
c     --- X - CN1 - CN5G - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN5G'
      tnlocal(isame,4) = 'wild'
c     --- X - CN2 - NN2G - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'wild'
c     --- X - CN2 - CN5 - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'wild'
c     --- X - CN5 - NN3A - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN3A'
      tnlocal(isame,4) = 'wild'
c     --- X - CN5 - NN3G - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN3G'
      tnlocal(isame,4) = 'wild'
c     --- X - CN5 - NN4 - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN4'
      tnlocal(isame,4) = 'wild'
c     --- X - CN5G - NN4 - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN5G'
      tnlocal(isame,3) = 'NN4'
      tnlocal(isame,4) = 'wild'
c     --- X - CN2 - NN3A - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN3A'
      tnlocal(isame,4) = 'wild'
c     --- X - CN2 - NN3A - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN3G'
      tnlocal(isame,4) = 'wild'
c     --- HN3B - CN3 - CN3B - HN3B  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN3B'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3B'
      tnlocal(isame,4) = 'HN3B'
c     --- CN1A - CN3 - CN3B - HN3B  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1A'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3B'
      tnlocal(isame,4) = 'HN3B'
c     --- CN1A - CN3 - CN3 - HN3  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1A'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- X - CN3 - CN3A - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3A'
      tnlocal(isame,4) = 'wild'
c     --- X - CN3 - CN3B - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3B'
      tnlocal(isame,4) = 'wild'
c     --- X - NN2 - CN3B - X  [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3B'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.15d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (8.50d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- X - CEL1 - CEL1 - X [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CEL1'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.900d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- X - CEL2 - CEL2 - X  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CEL2'
      tnlocal(isame,3) = 'CEL2'
      tnlocal(isame,4) = 'wild'
c     --- X - CE2 - CE2 - X  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CE2'
      tnlocal(isame,3) = 'CE2'
      tnlocal(isame,4) = 'wild'
c     -- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.200d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL2 - CEL1 - CEL2 - HEL2  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL2'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CEL2'
      tnlocal(isame,4) = 'HEL2'
c     --- CTL3 - CEL1 - CEL2 - HEL2  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL3'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CEL2'
      tnlocal(isame,4) = 'HEL2'
c     --- HEL1 - CEL1 - CEL2 - HEL2  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CEL2'
      tnlocal(isame,4) = 'HEL2'
c     --- CN5 - NN4 - CN4 - HN3  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'NN4'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'HN3'
c     --- CN5 - NN2 - CN4 - HN3  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'HN3'
c     --- CT3 - CE1 - CE2 - HE2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CE2'
      tnlocal(isame,4) = 'HE2'
c     --- HE1 - CE1 - CE2 - HE2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HE1'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CE2'
      tnlocal(isame,4) = 'HE2'
c     --- HE2 - CE2 - CE1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HE2'
      tnlocal(isame,2) = 'CE2'
      tnlocal(isame,3) = 'CE1'
      tnlocal(isame,4) = 'CT2'
c     --- X - CE1 - CE1 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CE1'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.030d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL1 - CEL1 - CTL2 - HAL2  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'HAL2'
c     --- CEL1 - CEL1 - CTL3 - HAL3  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL3'
      tnlocal(isame,4) = 'HAL3'
c     --- CE1 - CE1 - CT3 - HA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CE1'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.50d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (1.30d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL1 - CEL1 - CTL2 - CTL2  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CTL2'
c     --- CEL1 - CEL1 - CTL2 - CTL3  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CTL3'
c     --- CEL2 - CEL1 - CTL2 - CTL2  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL2'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CTL2'
c     --- CEL2 - CEL1 - CTL2 - CTL3  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL2'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CTL3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.120d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL2 - CEL1 - CTL2 - HAL2  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL2'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'HAL2'
c     --- HEL1 - CEL1 - CTL2 - CTL2  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CTL2'
c     --- HEL1 - CEL1 - CTL2 - CTL3  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CTL3'
c     --- CE2 - CE1 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CE2'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- HE1 - CE1 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HE1'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.050d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL2 - CEL1 - CTL3 - HAL3  [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL2'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL3'
      tnlocal(isame,4) = 'HAL3'
c     --- CE2 - CE1 - CT3 - HA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CE2'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (1.20d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.40d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      tclocal(7) = (1.30d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL2 - CEL1 - CTL2 - CEL1 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL2'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CEL1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(1) = (1.00d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.10d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.30d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      tclocal(10) = (0.20d0)*(CON_KCALTOK)
      tclocal(11) = 4.0d0
      tclocal(12) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CEL1 - CEL1 - CTL2 - CEL1 [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CEL1'
      tnlocal(isame,2) = 'CEL1'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CEL1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,5)
      tclocal(1) = (0.21d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (1.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.39d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.35d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      tclocal(10) = (0.11d0)*(CON_KCALTOK)
      tclocal(11) = 4.0d0
      tclocal(12) = (0.0d0)*(onepi/180.0d0)
      tclocal(13) = (0.09d0)*(CON_KCALTOK)
      tclocal(14) = 6.0d0
      tclocal(15) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CTL3 - CTL2 - CTL2 - CL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL3'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CL'
c     --- CTL2 - CTL2 - CTL2 - CL [par_all27_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CTL2'
      tnlocal(isame,2) = 'CTL2'
      tnlocal(isame,3) = 'CTL2'
      tnlocal(isame,4) = 'CL'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     *** Protein sections for dihedral terms (if not listed above)
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
      isame = 0
c     --- C - N - CP1 - C  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'C'
      tnlocal(isame,2) = 'N'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'C'
c     --- CC - CP1 - N - C  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CC'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'C'
c     --- CP2 - CP1 - N - C  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP2'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'C'
c     --- HB - CP1 - N - C  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HB'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'C'
c     --- CN7B - ON6 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'ON6'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8B'
c     --- HN5 - ON5 - CN7B - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN5'
      tnlocal(isame,2) = 'ON5'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'CN7B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CA - CA - CA - CA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- CA - CPT - CPT - CA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CA'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CA'
c     --- CPT - CA - CA - CA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CPT'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- CPT - CPT - CA - CA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CPT'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- CT2 - CA - CA - CA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- CT3 - CA - CA - CA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- OC - CA - CA - CA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OC'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- OH1 - CA - CA - CA  [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OH1'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)
      
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
      isame = 0
c     --- CA - CT2 - CT1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CA'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'C'
c     --- CC - CT1 - CT2 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CC'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CA'
c     --- HA - CT1 - CT2 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CA'
c     --- HA - CT3 - CT2 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CA'
c     --- N - CT1 - CT2 - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CA - CY - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CA'
      tnlocal(isame,2) = 'CY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CA'
c     --- CA - NY - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CA'
      tnlocal(isame,2) = 'NY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CA'
c     --- CT2 - CPH1 - NR1 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- CT2 - CPH1 - NR2 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR2'
      tnlocal(isame,4) = 'CPH2'
c     --- CT2 - CY - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CA'
c     --- CT2 - CY - CPT - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CPT'
c     --- CT3 - CPH1 - NR1 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- CY - CPT - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CY'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- H - NR3 - CPH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NR3'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT2'
c     --- H - NR3 - CPH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NR3'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT3'
c     --- HA - CA - CPT - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CPT'
c     --- HA - CA - NY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'NY'
      tnlocal(isame,4) = 'CPT'
c     --- HA - CY - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CA'
c     --- HA - CY - CPT - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CPT'
c     --- HP - CA - CA - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CPT'
c     --- HP - CA - CPT - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CPT'
c     --- HP - CA - CPT - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CY'
c     --- HR1 - CPH1 - NR3 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR3'
      tnlocal(isame,4) = 'H'
c     --- HR1 - CPH2 - NR1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR1'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH1'
c     --- HR1 - CPH2 - NR2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR1'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR2'
      tnlocal(isame,4) = 'CPH1'
c     --- HR2 - CPH2 - NR3 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR2'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR3'
      tnlocal(isame,4) = 'CPH1'
c     --- HR3 - CPH1 - NR1 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- HR3 - CPH1 - NR2 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR2'
      tnlocal(isame,4) = 'CPH2'
c     --- NR1 - CPH1 - CPH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT2'
c     --- NR1 - CPH1 - CPH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT3'
c     --- NR1 - CPH1 - CPH1 - HR3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'HR3'
c     --- NR2 - CPH1 - CPH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT2'
c     --- NR2 - CPH1 - CPH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT3'
c     --- NR2 - CPH1 - CPH1 - HR3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'HR3'
c     --- NY - CPT - CA - HP [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NY'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'HP'
c     --- X - CPB - C - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CPB'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'wild'
c     --- HN3 - CN3 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN3'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- CN1 - CN3T - CN3 - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1'
      tnlocal(isame,2) = 'CN3T'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'NN2B'
c     --- ON1C - CN1 - NN2 - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1C'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'HN2'
c     --- NN2 - CN3 - CN3D - CN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3D'
      tnlocal(isame,4) = 'CN2'
c     --- HN2 - NN2 - CN3B - HN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN2'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3B'
      tnlocal(isame,4) = 'HN3B'
c     --- CN3 - CN3A - CN3 - CN1A [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN3A'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'CN1A'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CP1 - C - N - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP1'
c     --- CP3 - N - C - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP3'
      tnlocal(isame,2) = 'N'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CP1'
c     --- CT1 - C - N - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP1'
c     --- CT1 - C - N - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP3'
c     --- CT2 - C - N - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP1'
c     --- CT2 - C - N - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP3'
c     --- CT3 - C - N - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP1'
c     --- CT3 - C - N - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP3'
c     --- O - C - N - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP3'
c     --- O - C - N - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'N'
      tnlocal(isame,4) = 'CP1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CP2 - CP3 - NP - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP2'
      tnlocal(isame,2) = 'CP3'
      tnlocal(isame,3) = 'NP'
      tnlocal(isame,4) = 'CP1'
c     --- CP3 - NP - CP1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP3'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'C'
c     --- CP3 - NP - CP1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP3'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CC'
c     --- CP3 - NP - CP1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP3'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CD'
c     --- CP3 - NP - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CP3'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CP2'
c     --- HA - CP3 - NP - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CP3'
      tnlocal(isame,3) = 'NP'
      tnlocal(isame,4) = 'CP1'
c     --- HB - CP1 - NP - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HB'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'NP'
      tnlocal(isame,4) = 'CP3'
c     --- HC - NP - CP1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'C'
c     --- HC - NP - CP1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CC'
c     --- HC - NP - CP1 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CD'
c     --- HC - NP - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CP2'
c     --- HC - NP - CP1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'HB'
c     --- HC - NP - CP3 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP3'
      tnlocal(isame,4) = 'CP2'
c     --- HC - NP - CP3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'NP'
      tnlocal(isame,3) = 'CP3'
      tnlocal(isame,4) = 'HA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CPH2 - NR1 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CPH2'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- CPH2 - NR2 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CPH2'
      tnlocal(isame,2) = 'NR2'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- NR1 - CPH2 - NR2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR2'
      tnlocal(isame,4) = 'CPH1'
c     --- NR2 - CPH1 - CPH1 - NR1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'NR1'
c     --- NR2 - CPH2 - NR1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR2'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH1'
c     --- CN5 - NN4 - CN4 - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'NN4'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'NN2'
c     --- NN2 - CN5 - CN5 - NN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'NN4'
c     --- ON1 - CN1 - CN5G - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN5G'
      tnlocal(isame,4) = 'CN5'
c     --- ON1 - CN1 - NN2G - CN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'CN2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CPH2 - NR3 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CPH2'
      tnlocal(isame,2) = 'NR3'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- NR3 - CPH1 - CPH1 - NR3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'NR3'
c     --- NR3 - CPH2 - NR3 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR3'
      tnlocal(isame,4) = 'CPH1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CPT - CPT - CY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CPT'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CA'
c     --- HA - CA - CPT - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CY'
c     --- NY - CA - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NY'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CPT'
c     --- NY - CPT - CA - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NY'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'HA'
c     --- X - C - C - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'wild'
c     --- NN2 - CN3 - CN3D - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3D'
      tnlocal(isame,4) = 'CN9'
c     --- CN2 - CN3D - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'CN3D'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- CN9 - CN3D - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN9'
      tnlocal(isame,2) = 'CN3D'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- CN4 - NN3A - CN2 - NN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN3A'
      tnlocal(isame,3) = 'CN2'
      tnlocal(isame,4) = 'NN1'
c     --- CN5 - CN5 - CN2 - NN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN2'
      tnlocal(isame,4) = 'NN1'
c     --- NN1 - CN2 - NN2G - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN1'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'CN1'
c     --- NN1 - CN2 - NN3G - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN1'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN3G'
      tnlocal(isame,4) = 'CN5'
c     --- CN3B - NN2 - CN3B - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3B'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3B'
      tnlocal(isame,4) = 'CN3'
c     --- HN2 - NN2 - CN3C - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN2'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3C'
      tnlocal(isame,4) = 'HN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CPT - CPT - NY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CPT'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'NY'
      tnlocal(isame,4) = 'CA'
c     --- CY - CA - NY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CY'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'NY'
      tnlocal(isame,4) = 'CPT'
c     --- NY - CPT - CPT - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NY'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CY'
c     --- CN1A - CN3 - CN3A - HN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1A'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3A'
      tnlocal(isame,4) = 'HN3B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CT1 - C - NH1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT1'
c     --- CT1 - NH1 - C - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CP1'
c     --- CT2 - C - NH1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT1'
c     --- CT2 - C - NH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT2'
c     --- CT2 - C - NH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT3'
c     --- CT2 - NH1 - C - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CP1'
c     --- CT2 - NH1 - C - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CT1'
c     --- CT3 - NH1 - C - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CT1'
c     --- CT3 - NH1 - C - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CT2'
c     --- CT3 - NH1 - C - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CT3'
c     --- CT3 - NH1 - C - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CP1'
c     --- CT3 - NH1 - C - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CT1'
c     --- CT3 - C - NH1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT1'
c     --- CT3 - C - NH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT2'
c     --- CT3 - C - NH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CT1 - CT1 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- CT2 - CT1 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- CT2 - CT2 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- CT3 - CT1 - NH1 - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'C'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CT1 - CT2 - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- CT1 - CT2 - CY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CA'
c     --- CT1 - CT2 - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CPT'
c     --- CT3 - CT2 - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- CT3 - CT2 - CY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CA'
c     --- CT3 - CT2 - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CPT'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
c     --- CT1 - CT2 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT1'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- CT3 - CT2 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
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
      isame = 0
c     --- CT2 - CPH1 - NR3 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR3'
      tnlocal(isame,4) = 'CPH2'
c     --- H - NH1 - C - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CP1'
c     --- H - NH1 - C - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CT1'
c     --- H - NH1 - C - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CT2'
c     --- H - NH1 - C - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'CT3'
c     --- H - NH2 - CC - CP1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH2'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'CP1'
c     --- HA - CA - CA - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'HA'
c     --- HR1 - CPH1 - NR3 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR3'
      tnlocal(isame,4) = 'CPH2'
c     --- NR3 - CPH1 - CPH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT2'
c     --- NR3 - CPH1 - CPH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT3'
c     --- NR3 - CPH1 - CPH1 - HR1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'HR1'
c     --- O - C - NH1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT1'
c     --- O - C - NH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT2'
c     --- O - C - NH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'CT3'
c     --- O - C - NH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'H'
c     --- X - NN1 - CN1A - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'NN1'
      tnlocal(isame,3) = 'CN1A'
      tnlocal(isame,4) = 'wild'
c     --- CN1A - CN3 - CN3C - NN2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1A'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3C'
      tnlocal(isame,4) = 'NN2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CT2 - CT2 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- NH1 - CT1 - C - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'N'
c     --- NH1 - CT2 - C - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'N'
c     --- NH3 - CT1 - C - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH3'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'N'
c     --- NH3 - CT1 - CC - NH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH3'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'NH2'
c     --- NH3 - CT2 - C - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'N'
c     --- NH3 - CT2 - C - NH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'NH1'
c     --- NH3 - CT2 - CC - NH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'NH2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CT2 - CT2 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- CT3 - CT2 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- CT3 - CT2 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT3'
c     --- CN9 - CN8 - CN8 - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN9'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN9'
c     --- CN9 - CN8 - CN8 - CN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN9'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN8'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CT2 - SM - SM - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'SM'
      tnlocal(isame,3) = 'SM'
      tnlocal(isame,4) = 'CT2'
c     --- CT3 - SM - SM - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'SM'
      tnlocal(isame,3) = 'SM'
      tnlocal(isame,4) = 'CT3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CT3 - CT2 - S - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'CT3'
c     --- CT3 - S - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (10.00d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CY - CPT - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CY'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CA'
c     --- NY - CPT - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NY'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CA'
c     --- CN3D - CN2 - NN3 - CN1 [par_all27_na.inp]
      tnlocal(isame,1) = 'CN3D'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN3'
      tnlocal(isame,4) = 'CN1'
c     --- CN5 - CN2 - NN3A - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN3A'
      tnlocal(isame,4) = 'CN4'
c     --- NN2B - CN5 - CN5G - NN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2B'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN5G'
      tnlocal(isame,4) = 'NN4'
c     --- NN3G - CN5 - CN5G - NN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3G'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN5G'
      tnlocal(isame,4) = 'NN4'
c     --- CN1 - CN5G - CN5 - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1'
      tnlocal(isame,2) = 'CN5G'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'NN2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- H - NH2 - CC - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH2'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'CT1'
c     --- H - NH2 - CC - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH2'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'CT2'
c     --- H - NH2 - CC - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NH2'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'CT3'
c     --- H - NR3 - CPH1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NR3'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CPH1'
c     --- HA - CC - NH2 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'NH2'
      tnlocal(isame,4) = 'H'
c     --- NR3 - CPH2 - NR3 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NR3'
      tnlocal(isame,2) = 'CPH2'
      tnlocal(isame,3) = 'NR3'
      tnlocal(isame,4) = 'H'
c     --- O - CC - NH2 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'NH2'
      tnlocal(isame,4) = 'H'
c     --- HA - C- NH1 - H (formamide cap on gramicidin)
c     --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c         `High-resolution polypeptide structure in a lamellar 
c          phase lipid environment from solid state NMR derived 
c          orientational constraints'
c         Structure, 5, 1655-1669 (1997).
c     --- T.B. Woolf personal communication 08-16-2001 (M.G. Martin)
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NH1'
      tnlocal(isame,4) = 'H'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- H - NY - CA - CY [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NY'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CY'
c     --- H - NY - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CA'
c     --- H - NY - CPT - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'NY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CPT'
c     --- X - CN2 - CN3 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- H - OH1 - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'OH1'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- H - OH1 - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'OH1'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT1'
c     --- H - OH1 - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'OH1'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT3'
c     --- CN7 - CN8B - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN8B'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- H - OH1 - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'OH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT1'
c     --- H - OH1 - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'OH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- H - OH1 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'H'
      tnlocal(isame,2) = 'OH1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HA - CA - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- HA - CA - CA - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CPT'
c     --- NY - CA - CY - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NY'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CT2'
c     --- NY - CA - CY - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NY'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'HA'
c     --- NY - CA - CY - HP [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NY'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'HP'
c     --- X - CN4 - NN3I - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN3I'
      tnlocal(isame,4) = 'wild'
c     --- X - CN4 - NN3A - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN3A'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HA - CA - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CPT'
c     --- HA - CA - CY - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CT2'
c     --- HA - CY - CA - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CY'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CPT'
c     --- HA - CY - CA - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CY'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'HA'
c     --- HP - CA - CY - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CT2'
c     --- CN5 - CN5 - NN2 - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'HN2'
c     --- NN4 - CN4 - NN2 - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN4'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'HN2'
c     --- NN2G - CN2 - NN1 - HN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2G'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN1'
      tnlocal(isame,4) = 'HN1'
c     --- NN3G - CN2 - NN1 - HN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3G'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN1'
      tnlocal(isame,4) = 'HN1'
c     --- HN2 - NN2B - CN5 - CN5G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN2'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'CN5G'
c     --- HN2 - NN2B - CN5 - NN3G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN2'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'NN3G'
c     --- HN2 - NN2B - CN4 - NN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN2'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'NN4'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HA - CT2 - CY - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CA'
c     --- HA - CT2 - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CPT'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HA - CT2 - S - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'CT3'
c     --- HA - CT3 - S - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'S'
      tnlocal(isame,4) = 'CT2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HC - NH2 - CT2 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'NH2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HB'
c     --- HC - NH2 - CT2 - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'NH2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CD'
c     --- X - CT3 - NH2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'NH2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HP - CA - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- HP - CA - CA - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CT2'
c     --- HP - CA - CA - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CT3'
c     --- OC - CA - CA - HP [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OC'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'HP'
c     --- OH1 - CA - CA - HP [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OH1'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'HP'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HP - CA - CA - HP [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'HP'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HP - CA - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'CY'
      tnlocal(isame,4) = 'CPT'
c     --- HP - CA - NY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'NY'
      tnlocal(isame,4) = 'CPT'
c     --- HR3 - CPH1 - CPH1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT2'
c     --- HR3 - CPH1 - CPH1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'CT3'
c     --- HR3 - CPH1 - CPH1 - HR3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'HR3'
c     --- NN1 - CN2 - NN3 - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN1'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN3'
      tnlocal(isame,4) = 'CN1'
c     --- NN1 - CN2 - CN3 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN1'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'CN3'
c     --- NN1 - CN2 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN1'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- X - CN2 - NN3 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN3'
      tnlocal(isame,4) = 'wild'
c     --- CN3D - CN2 - NN1 - HN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3D'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN1'
      tnlocal(isame,4) = 'HN1'
c     --- NN3A - CN4 - NN3A - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3A'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN3A'
      tnlocal(isame,4) = 'CN5'
c     --- CN4 - NN3A - CN5 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN3A'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'CN5'
c     --- NN3A - CN2 - CN5 - NN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3A'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'NN4'
c     --- CN2 - CN5 - NN4 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN4'
      tnlocal(isame,4) = 'CN4'
c     --- CN4 - NN3A - CN5 - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN3A'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'NN2'
c     --- NN2G - CN2 - NN3G - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2G'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN3G'
      tnlocal(isame,4) = 'CN5'
c     --- NN3G - CN5 - CN5G - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3G'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN5G'
      tnlocal(isame,4) = 'CN1'
c     --- NN2G - CN1 - CN5G - NN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2G'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN5G'
      tnlocal(isame,4) = 'NN4'
c     --- CN1 - CN5G - NN4 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1'
      tnlocal(isame,2) = 'CN5G'
      tnlocal(isame,3) = 'NN4'
      tnlocal(isame,4) = 'CN4'
c     --- CN2 - NN3G - CN5 - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'NN3G'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'NN2B'
c     --- NN3G - CN5 - NN2B - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3G'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN4'
c     --- X - CN4 - NN4 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN4'
      tnlocal(isame,4) = 'wild'
c     --- HN3B - CN3 - CN3A - HN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN3B'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3A'
      tnlocal(isame,4) = 'HN3B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HP - CA - NY - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CA'
      tnlocal(isame,3) = 'NY'
      tnlocal(isame,4) = 'H'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HP - CY - CPT - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CA'
c     --- HP - CY - CPT - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HP'
      tnlocal(isame,2) = 'CY'
      tnlocal(isame,3) = 'CPT'
      tnlocal(isame,4) = 'CPT'
c     --- NY - CPT - CA - CA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NY'
      tnlocal(isame,2) = 'CPT'
      tnlocal(isame,3) = 'CA'
      tnlocal(isame,4) = 'CA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HS - S - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HS'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT1'
c     --- HS - S - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HS'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- HS - S - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HS'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- HS - S - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HS'
      tnlocal(isame,2) = 'S'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- X - CT1 - CT1 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'wild'
c     --- X - CT1 - CT2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'wild'
c     --- X - CT1 - CT3 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- N - C - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CP2'
c     --- NH1 - C - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CP2'
c     --- NH2 - CC - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH2'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CP2'
c     --- O - C - CP1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'HB'
c     --- O - CC - CP1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'HB'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- N - C - CP1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'HB'
c     --- O - C - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CP2'
c     --- O - CC - CP1 - CP2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'CP2'
c     --- NH1 - C - CP1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'HB'
c     --- NH2 - CC - CP1 - HB [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH2'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'HB'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- N - C - CP1 - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'N'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'N'
c     --- NH1 - C - CP1 - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'N'
c     --- NH2 - CC - CP1 - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH2'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'N'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- NH1 - C - CT2 - NH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'NH1'
c     --- NH3 - CT1 - C - NH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH3'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'NH1'
c     --- NH1 - C - CT1 - NH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NH1'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'NH1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- NP - CP1 - C - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NP'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'N'
c     --- NP - CP1 - C - NH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NP'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'NH1'
c     --- NP - CP1 - CC - NH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NP'
      tnlocal(isame,2) = 'CP1'
      tnlocal(isame,3) = 'CC'
      tnlocal(isame,4) = 'NH2'
c     --- ON6B - CN7B - NR1 - CPH2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- ON6 - CN7B - NR1 - CPH2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- O - C - CP1 - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'N'
c     --- O - CC - CP1 - N [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CP1'
      tnlocal(isame,4) = 'N'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- O - C - CT1 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT1'
c     --- O - C - CT1 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT2'
c     --- O - C - CT1 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'CT3'
c     --- O - C - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'O'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- ON6 - CN7B - NN2B - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN4'
c     --- ON6B - CN7B - NN2B - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN4'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- OC - CC - CT1 - NH3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OC'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CT1'
      tnlocal(isame,4) = 'NH3'
c     --- OC - CC - CT2 - NH3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'OC'
      tnlocal(isame,2) = 'CC'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'NH3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- S - CT2 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'S'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- SM - CT2 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'SM'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- SM - SM - CT2 - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'SM'
      tnlocal(isame,2) = 'SM'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT1'
c     --- SM - SM - CT2 - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'SM'
      tnlocal(isame,2) = 'SM'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- SM - SM - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'SM'
      tnlocal(isame,2) = 'SM'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- SM - SM - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'SM'
      tnlocal(isame,2) = 'SM'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- SS - CS - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'SS'
      tnlocal(isame,2) = 'CS'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- X - CN8 - CN8 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'wild'
c     --- X - CN8 - CN9 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - C - NC2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'NC2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CT1 - OS - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT1'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'wild'
c     --- X - CT2 - OS - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'wild'
c     --- X - CT3 - OS - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'OS'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CT2 - CT2 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'wild'
c     --- HN7 - CN7 - CN7B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON5'
c     --- HN7 - CN7 - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- HN7 - CN7 - CN7 - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7B'
c     --- HN7 - CN7B - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- HN7 - CN7B - NN2 - CN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3B'
c     --- HN7 - CN7B - NN2 - CN3C [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3C'
c     --- HN7 - CN7 - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- HN7 - CN7 - CN7B - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6'
c     --- HN7 - CN7 - CN7B - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON2'
c     --- HN7 - CN7B - CN7B - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON2'
c     --- HN7 - CN7 - CN7B - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NN2'
c     --- HN7 - CN7B - NR1 - CPH2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- HN7 - CN7B - NN2 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN4'
c     --- HN7 - CN7B - NN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3'
c     --- HN7 - CN7B - NN2B - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN3'
c     --- HN7 - CN7B - NN2B - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN4'
c     --- HN8 - CN8B - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- HN7 - CN7 - CN8B - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'ON2'
c     --- HN7 - CN7 - CN8B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'ON5'
c     --- ON2 - CN7 - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'HN8'
c     --- ON2 - CN7 - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- ON2 - CN7 - CN7C - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'HN7'
c     --- HN7 - CN7 - CN7 - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON2'
c     --- HN7 - CN7 - CN7 - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON5'
c     --- HN7 - CN7 - CN8B - HN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'HN8'
c     --- HN7 - CN7 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8B'
c     --- HN7 - CN7 - CN8 - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN7B'
c     --- HN7 - CN7B - CN8 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN7'
c     --- HN8 - CN8 - CN7B - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6'
c     --- HN7 - CN7 - CN7 - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'HN7'
c     --- HN7 - CN7B - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'HN8'
c     --- HN7 - CN7 - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'HN8'
c     --- HN8 - CN8 - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- HN7 - CN7 - ON6 - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON6'
      tnlocal(isame,4) = 'CN7B'
c     --- HN8 - CN8 - CN7B - ON6B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6B'
c     --- HN7 - CN7 - ON6B - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON6B'
      tnlocal(isame,4) = 'CN7B'
c     --- HN7 - CN7 - CN7B - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'CN7B'
c     --- HN7 - CN7B - CN7B - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'CN7'
c     --- HN7 - CN7B - CN7B - ON6B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6B'
c     --- HN7 - CN7 - CN7C - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'CN7B'
c     --- HN7 - CN7B - CN7C - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'CN7'
c     --- HN7 - CN7C - CN7B - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6'
c     --- HN7 - CN7B - CN7C - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'HN7'
c     --- HN7 - CN7 - CN7C - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'HN7'
c     --- HN7 - CN7C - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CT3 - CT3 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- X - CT2 - NH3 - X [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'NH3'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     *** Nucleic acid dihedral terms (if not listed above) ***
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(1) = (1.20d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.10d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.10d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      tclocal(10) = (0.00d0)*(CON_KCALTOK)
      tclocal(11) = 6.0d0
      tclocal(12) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - P - ON2 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN7'
c     --- ON2 - P - ON2 - CN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN8'
c     --- ON2 - P - ON2 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN8B'
c     --- ON2 - P - ON2 - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'P'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'CN9'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.500d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- NN2B - CN1T - NN2U - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2B'
      tnlocal(isame,2) = 'CN1T'
      tnlocal(isame,3) = 'NN2U'
      tnlocal(isame,4) = 'CN1'
c     --- CN1T - NN2U - CN1 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1T'
      tnlocal(isame,2) = 'NN2U'
      tnlocal(isame,3) = 'CN1'
      tnlocal(isame,4) = 'CN3'
c     --- NN2U - CN1 - CN3 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2U'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'CN3'
c     --- CN3 - CN3 - NN2B - CN1T [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN1T'
c     --- CN3 - NN2B - CN1T - NN2U [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN1T'
      tnlocal(isame,4) = 'NN2U'
c     --- HN2 - NN2B - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN2'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- NN3I - CN4 - NN2G - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3I'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'HN2'
c     --- NN3I - CN4 - NN2G - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3I'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'CN1'
c     --- HN3 - CN4 - NN2G - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN3'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'HN2'
c     --- HN3 - CN4 - NN2G - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN3'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'CN1'
c     --- X - CN4 - NN2B - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'wild'
c     --- X - CN5 - NN2 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'wild'
c     --- X - CN5 - NN2B - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'wild'
c     --- X - CN4 - NN2 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (6.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - CN3 - CN3 - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'NN2B'
c     --- HN3 - CN3 - CN1 - ON1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN3'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN1'
      tnlocal(isame,4) = 'ON1'
c     --- CN1 - NN3 - CN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1'
      tnlocal(isame,2) = 'NN3'
      tnlocal(isame,3) = 'CN2'
      tnlocal(isame,4) = 'CN3'
c     --- CN2 - CN3 - CN3 - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'NN2'
c     --- CN5 - CN5 - NN4 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN4'
      tnlocal(isame,4) = 'CN4'
c     --- NN4 - CN4 - NN2 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN4'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN5'
c     --- CN4 - NN2 - CN5 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'CN5'
c     --- CN5 - CN5G - NN4 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'CN5G'
      tnlocal(isame,3) = 'NN4'
      tnlocal(isame,4) = 'CN4'
c     --- NN4 - CN4 - NN2B - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN4'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN5'
c     --- CN4 - NN2B - CN5 - CN5G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'CN5G'
c     --- CN5 - CN5 - NN3I - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'NN3I'
      tnlocal(isame,4) = 'CN4'
c     --- NN2G - CN4 - NN3I - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2G'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN3I'
      tnlocal(isame,4) = 'CN5'
c     --- CN3 - CN3A - CN3 - CN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN3A'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'CN3B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.800d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- NN2B - CN1T - NN2U - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2B'
      tnlocal(isame,2) = 'CN1T'
      tnlocal(isame,3) = 'NN2U'
      tnlocal(isame,4) = 'HN2'
c     --- CN3 - CN1 - NN2U - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2U'
      tnlocal(isame,4) = 'HN2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.600d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - CN3 - NN2B - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'HN2'
c     --- NN2U - CN1T - NN2B - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2U'
      tnlocal(isame,2) = 'CN1T'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'HN2'
c     --- ON1C - CN1 - NN3 - CN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1C'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN3'
      tnlocal(isame,4) = 'CN2'
c     --- ON1C - CN1 - NN3 - CN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON1C'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.800d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1T - NN2B - CN3 - CN3T [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1T'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'CN3T'
c     --- NN2U - CN1 - CN3T - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2U'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN3T'
      tnlocal(isame,4) = 'CN3'
c     --- CN2 - NN3A - CN4 - NN3A [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'NN3A'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'NN3A'
c     --- CN4 - NN3A - CN5 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN4'
      tnlocal(isame,2) = 'NN3A'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'CN5'
c     --- CN5 - CN5 - CN2 - NN3A [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN2'
      tnlocal(isame,4) = 'NN3A'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.600d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- NN2B - CN3 - CN3T - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2B'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3T'
      tnlocal(isame,4) = 'CN9'
c     --- NN2B - CN1 - CN3T - CN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2B'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN3T'
      tnlocal(isame,4) = 'CN9'
c     --- HN3 - CN4 - NN4 - CN5G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN3'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN4'
      tnlocal(isame,4) = 'CN5G'
c     --- HN3 - CN4 - NN2B - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN3'
      tnlocal(isame,2) = 'CN4'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.460d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - CN3T - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1'
      tnlocal(isame,2) = 'CN3T'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- CN3 - CN3T - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN3T'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.800d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3T - CN1 - NN2U - HN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3T'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2U'
      tnlocal(isame,4) = 'HN2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.600d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - NN2 - CN1 - NN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN1'
      tnlocal(isame,4) = 'NN3'
c     --- NN2 - CN1 - NN3 - CN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN3'
      tnlocal(isame,4) = 'CN2'
c     --- NN3 - CN2 - CN3 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'CN3'
c     --- CN3 - CN3 - NN2 - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.400d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- NN3 - CN2 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- NN2 - CN3 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.600d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - CN3 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- CN3 - NN2 - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.350d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - CN3D - CN9 - HN9 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'CN3D'
      tnlocal(isame,3) = 'CN9'
      tnlocal(isame,4) = 'HN9'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (8.500d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - NN3A - CN4 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'NN3A'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'HN3'
c     --- CN5 - NN3A - CN4 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'NN3A'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'HN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5 - CN2 - NN1 - HN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN1'
      tnlocal(isame,4) = 'HN1'
c     --- NN3A - CN2 - NN1 - HN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3A'
      tnlocal(isame,2) = 'CN2'
      tnlocal(isame,3) = 'NN1'
      tnlocal(isame,4) = 'HN1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (7.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- NN3A - CN5 - CN5 - NN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN3A'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'NN4'
c     --- CN2 - CN5 - CN5 - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'CN5'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'NN2'
c     --- CN3 - NN2 - CN3B - HN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3B'
      tnlocal(isame,4) = 'HN3B'
c     --- NN2 - CN3B - CN3 - HN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2'
      tnlocal(isame,2) = 'CN3B'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3B'
c     --- CN3 - CN1A - CN3 - HN3B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN1A'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3B'
c     --- CN3 - CN1A - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'CN1A'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- NN2 - CN3B - CN3 - CN3A [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2'
      tnlocal(isame,2) = 'CN3B'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'CN3A'
c     --- CN3 - NN2 - CN3C - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3C'
      tnlocal(isame,4) = 'HN3'
c     --- NN2 - CN3C - CN3 - HN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'NN2'
      tnlocal(isame,2) = 'CN3C'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      isame = 0
c     --- CN1 - NN2G - CN2 - NN3G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN1'
      tnlocal(isame,2) = 'NN2G'
      tnlocal(isame,3) = 'CN2'
      tnlocal(isame,4) = 'NN3G'
c     --- CN2 - NN3G - CN5 - CN5G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN2'
      tnlocal(isame,2) = 'NN3G'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'CN5G'
c     --- CN5 - CN5G - CN1 - NN2G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5'
      tnlocal(isame,2) = 'CN5G'
      tnlocal(isame,3) = 'CN1'
      tnlocal(isame,4) = 'NN2G'
c     --- CN5G - CN1 - NN2G - CN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5G'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'CN2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (16.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN5G - NN4 - CN4 - NN2B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN5G'
      tnlocal(isame,2) = 'NN4'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'NN2B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.600d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HN2 - NN2G - CN1 - CN5G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN2'
      tnlocal(isame,2) = 'NN2G'
      tnlocal(isame,3) = 'CN1'
      tnlocal(isame,4) = 'CN5G'
c     --- HN2 - NN2G - CN2 - NN3G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN2'
      tnlocal(isame,2) = 'NN2G'
      tnlocal(isame,3) = 'CN2'
      tnlocal(isame,4) = 'NN3G'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.900d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- X - CN1 - NN2 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'wild'
c     --- X - CN1T - NN2B - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1T'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'wild'
c     --- X - CN1 - NN2G - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2G'
      tnlocal(isame,4) = 'wild'
c     --- X - CN1 - NN2U - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'NN2U'
      tnlocal(isame,4) = 'wild'
c     --- X - CN1T - NN2U - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1T'
      tnlocal(isame,3) = 'NN2U'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.380d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3A - CN3 - CN1A - ON1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3A'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN1A'
      tnlocal(isame,4) = 'ON1'
c     --- CN3B - CN3 - CN1A - ON1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3B'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN1A'
      tnlocal(isame,4) = 'ON1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.350d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.620d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3A - CN3 - CN1A - NN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3A'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN1A'
      tnlocal(isame,4) = 'NN1'
c     --- CN3B - CN3 - CN1A - NN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3B'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN1A'
      tnlocal(isame,4) = 'NN1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.000d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- X - CN3 - CN8 - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (1.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (1.000d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.400d0)*(CON_KCALTOK)
      tclocal(8) = 6.0d0
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8 - CN3 - CN1A - ON1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN1A'
      tnlocal(isame,4) = 'ON1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.350d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.400d0)*(CON_KCALTOK)
      tclocal(8) = 6.0d0
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8 - CN3 - CN1A - NN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN1A'
      tnlocal(isame,4) = 'NN1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.300d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (1.950d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3C - CN3 - CN1A - ON1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3C'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN1A'
      tnlocal(isame,4) = 'ON1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (1.100d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (1.950d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3C - CN3 - CN1A - NN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3C'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN1A'
      tnlocal(isame,4) = 'NN1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8 - CN3 - CN3C - NN2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3C'
      tnlocal(isame,4) = 'NN2'
c     --- X - CN3 - CN3C - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN3C'
      tnlocal(isame,4) = 'wild'
c     --- X - NN2 - CN3C - X [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3C'
      tnlocal(isame,4) = 'wild'
c     --- CN3C - NN2 - CN3C - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3C'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3C'
      tnlocal(isame,4) = 'CN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.000d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3C - CN3 - CN8 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN3C'
      tnlocal(isame,2) = 'CN3'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.030d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.030d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- P2 - ON2 - P - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P2'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'P'
      tnlocal(isame,4) = 'ON2'
c     --- P - ON2 - P2 - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'P2'
      tnlocal(isame,4) = 'ON2'
c     --- P2 - ON2 - P2 - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P2'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'P2'
      tnlocal(isame,4) = 'ON2'
c     --- P - ON2 - P - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'P'
      tnlocal(isame,4) = 'ON2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.030d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- P - ON2 - P - ON3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'P'
      tnlocal(isame,4) = 'ON3'
c     --- P - ON2 - P2 - ON3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'P2'
      tnlocal(isame,4) = 'ON3'
c     --- P2 - ON2 - P - ON3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P2'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'P'
      tnlocal(isame,4) = 'ON3'
c     --- P2 - ON2 - P2 - ON3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'P2'
      tnlocal(isame,2) = 'ON2'
      tnlocal(isame,3) = 'P2'
      tnlocal(isame,4) = 'ON3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.400d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON6 - CN7B - CN7 - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON5'
c     --- ON6 - CN7 - CN8B - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'ON2'
c     --- ON6B - CN7 - CN8B - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'ON2'
c     --- ON6 - CN7 - CN8B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'ON5'
c     --- ON6B - CN7 - CN8B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'ON5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.600d0)*(CON_KCALTOK)
      tclocal(2) = 6.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN7B - ON6 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'ON6'
      tnlocal(isame,4) = 'CN7'
c     --- CN8 - CN7B - ON6 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'ON6'
      tnlocal(isame,4) = 'CN7'
c     --- CN7C - CN7B - ON6 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'ON6'
      tnlocal(isame,4) = 'CN7'
c     --- CN7 - CN7 - ON6 - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON6'
      tnlocal(isame,4) = 'CN7B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.400d0)*(CON_KCALTOK)
      tclocal(2) = 6.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - CN7 - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- CN7 - CN7B - CN7B - ON6B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6B'
c     --- CN7B - CN8 - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- CN7B - CN7C - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.600d0)*(CON_KCALTOK)
      tclocal(2) = 6.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN7 - CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6'
c     --- CN7 - CN7 - CN7B - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6'
c     --- CN7 - CN8 - CN7B - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6'
c     --- CN7 - CN7C - CN7B - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7C'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.800d0)*(CON_KCALTOK)
      tclocal(2) = 6.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.400d0)*(CON_KCALTOK)
      tclocal(5) = 5.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (2.000d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - CN7 - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- ON2 - CN7B - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- ON5 - CN7 - CN7 - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN7'
c     --- ON2 - CN7 - CN8 - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN7B'
c     --- ON2 - CN7 - CN7C - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'CN7B'
c     --- ON5 - CN7 - CN7C - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'CN7B'
c     --- ON5 - CN7 - CN8 - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'CN7B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(2) = 6.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.300d0)*(CON_KCALTOK)
      tclocal(5) = 5.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.600d0)*(CON_KCALTOK)
      tclocal(8) = 4.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      tclocal(10) = (0.200d0)*(CON_KCALTOK)
      tclocal(11) = 3.0d0
      tclocal(12) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - CN7 - CN7B - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6'
c     --- ON2 - CN7 - CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6'
c     --- ON5 - CN7 - CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN7 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8B'
c     --- CN8 - CN7 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8B'
c     --- CN7C - CN7 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.200d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.500d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (-70.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.600d0)*(CON_KCALTOK)
      tclocal(8) = 1.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON6 - CN7B - NR1 - CPH1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH1'
c     --- ON6B - CN7B - NR1 - CPH1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.200d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8 - CN7B - NR1 - CPH1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH1'
c     --- CN7B - CN7B - NR1 - CPH1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.100d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.300d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.000d0)*(CON_KCALTOK)
      tclocal(8) = 1.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8 - CN7B - NR1 - CPH2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- CN7B - CN7B - NR1 - CPH2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.300d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - NN2 - CN4 - HN3[par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'HN3'
c     --- CN7B - NN2 - CN3 - HN3[par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- CN7B - NN2B - CN3 - HN3[par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'HN3'
c     --- CN7B - NN2B - CN4 - HN3[par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'HN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (11.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - NN2 - CN5 - CN5[par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'CN5'
c     --- CN7B - NN2 - CN4 - NN4[par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'NN4'
c     --- CN7B - NN2 - CN4 - NN3A [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'NN3A'
c     --- CN7B - NN2 - CN1 - ON1C [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN1'
      tnlocal(isame,4) = 'ON1C'
c     --- CN7B - NN2 - CN1 - NN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN1'
      tnlocal(isame,4) = 'NN3'
c     --- CN7B - NN2 - CN3 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2'
      tnlocal(isame,3) = 'CN3'
      tnlocal(isame,4) = 'CN3'
c     --- CN7B - NN2B - CN1T - ON1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN1T'
      tnlocal(isame,4) = 'ON1'
c     --- CN7B - NN2B - CN1T - NN2U [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN1T'
      tnlocal(isame,4) = 'NN2U'
c     --- CN7B - NN2B - CN1T - CN3T [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN1T'
      tnlocal(isame,4) = 'CN3T'
c     --- CN7B - NN2B - CN4 - NN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN4'
      tnlocal(isame,4) = 'NN4'
c     --- CN7B - NN2B - CN5 - CN5G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'CN5G'
c     --- CN7B - NN2B - CN5 - NN3G [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NN2B'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'NN3G'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.1000d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON6 - CN7B - NN2 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN5'
c     --- ON6B - CN7B - NN2 - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.1000d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON6 - CN7B - NN2 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN4'
c     --- ON6B - CN7B - NN2 - CN4 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN4'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.0000d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON6 - CN7B - NN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3'
c     --- ON6B - CN7B - NN2 - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN3'
c     --- CN7B - CN7B - NN2 - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN1'
c     --- CN8 - CN7B - NN2 - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN1'
c     --- CN7C - CN7B - NN2 - CN1 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2'
      tnlocal(isame,4) = 'CN1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.7000d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.8000d0)*(CON_KCALTOK)
      tclocal(5) = 1.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON6 - CN7B - NN2B - CN1T [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN1T'
c     --- ON6B - CN7B - NN2B - CN1T [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN1T'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.9000d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON6 - CN7B - NN2B - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN3'
c     --- ON6B - CN7B - NN2B - CN3 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.2000d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8 - CN7B - NN2B - CN1T [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN1T'
c     --- CN7B - CN7B - NN2B - CN1T [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN1T'
c     --- CN7C - CN7B - NN2B - CN1T [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN1T'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.2000d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (1.1000d0)*(CON_KCALTOK)
      tclocal(5) = 1.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON6 - CN7B - NN2B - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN5'
c     --- ON6B - CN7B - NN2B - CN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON6B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'NN2B'
      tnlocal(isame,4) = 'CN5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.2000d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (120.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN8B - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN8B'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'P'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(1) = (0.2000d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.8000d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.4000d0)*(CON_KCALTOK)
      tclocal(8) = 2.0d0
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      tclocal(10) = (2.5000d0)*(CON_KCALTOK)
      tclocal(11) = 1.0d0
      tclocal(12) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN7 - CN8B - ON2 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'ON2'
c     --- CN7 - CN7 - CN8B - ON5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8B'
      tnlocal(isame,4) = 'ON5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.195d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HN8 - CN8B - CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6'
c     --- HN8 - CN8B - CN7 - ON6B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6B'
c     --- HN8 - CN8 - CN8 - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN8'
      tnlocal(isame,2) = 'CN8'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'ON6'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.200d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - CN7 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.200d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.800d0)*(CON_KCALTOK)
      tclocal(5) = 3.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - CN7 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8B'
c     --- ON5 - CN7 - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(1) = (0.400d0)*(CON_KCALTOK)
      tclocal(2) = 6.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.000d0)*(CON_KCALTOK)
      tclocal(5) = 5.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.000d0)*(CON_KCALTOK)
      tclocal(8) = 4.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      tclocal(10) = (1.600d0)*(CON_KCALTOK)
      tclocal(11) = 3.0d0
      tclocal(12) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - CN7 - CN7 - ON8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6B'
c     --- ON2 - CN7B - CN7B - ON8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'ON6B'
c     --- ON5 - CN7 - CN7 - ON8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.000d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - ON6B - CN7 - CN8B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'ON6B'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'CN8B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.600d0)*(CON_KCALTOK)
      tclocal(2) = 6.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.000d0)*(CON_KCALTOK)
      tclocal(5) = 5.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (1.600d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- ON2 - CN7 - CN7B - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON2'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'CN7B'
c     --- ON5 - CN7 - CN7B - CN7B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'CN7B'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,5)
      tclocal(1) = (0.600d0)*(CON_KCALTOK)
      tclocal(2) = 5.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.200d0)*(CON_KCALTOK)
      tclocal(5) = 4.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      tclocal(7) = (0.000d0)*(CON_KCALTOK)
      tclocal(8) = 3.0d0
      tclocal(9) = (180.0d0)*(onepi/180.0d0)
      tclocal(10) = (0.400d0)*(CON_KCALTOK)
      tclocal(11) = 2.0d0
      tclocal(12) = (0.0d0)*(onepi/180.0d0)
      tclocal(13) = (1.900d0)*(CON_KCALTOK)
      tclocal(14) = 1.0d0
      tclocal(15) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN7 - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'P'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.500d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN8 - CN7 - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'P'
c     --- CN7B - CN7 - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'P'
c     --- CN7B - CN7B - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'P'
c     --- CN7 - CN7B - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7B'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'P'
c     --- CN7 - CN7B - ON2 - P [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON2'
      tnlocal(isame,4) = 'P'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.300d0)*(CON_KCALTOK)
      tclocal(5) = 2.0d0
      tclocal(6) = (180.0d0)*(onepi/180.0d0)
      tclocal(7) = (1.500d0)*(CON_KCALTOK)
      tclocal(8) = 1.0d0
      tclocal(9) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7 - CN7 - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- CN8 - CN7 - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN8'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.800d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.500d0)*(CON_KCALTOK)
      tclocal(5) = 1.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - CN7 - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- CN7C - CN7 - ON5 - HN5 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CN7C'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'ON5'
      tnlocal(isame,4) = 'HN5'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.195d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HN7 - CN7 - CN7 - ON6 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6'
c     --- HN7 - CN7 - CN7 - ON6B [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN7'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7'
      tnlocal(isame,4) = 'ON6B'
c     --- ON5 - CN7 - CN8 - HN8 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN8'
      tnlocal(isame,4) = 'HN8'
c     --- ON5 - CN7 - CN7B - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- ON5 - CN7 - CN7C - HN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'ON5'
      tnlocal(isame,2) = 'CN7'
      tnlocal(isame,3) = 'CN7C'
      tnlocal(isame,4) = 'HN7'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,2)
      tclocal(1) = (0.300d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(4) = (0.000d0)*(CON_KCALTOK)
      tclocal(5) = 1.0d0
      tclocal(6) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HN5 - CON5 - CN7B - CN7 [par_all27_na.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HN5'
      tnlocal(isame,2) = 'ON5'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'CN7'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0500d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- wild - FE - CM - wild [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'FE'
      tnlocal(isame,3) = 'CM'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CE2 - CE1 - CT2 - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'CE2'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.8700d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HE1 - CE1 - CT2 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HE1'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.3400d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HE1 - CE1 - CT3 - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      tnlocal(isame,1) = 'HE1'
      tnlocal(isame,2) = 'CE1'
      tnlocal(isame,3) = 'CT3'
      tnlocal(isame,4) = 'HA'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffcharmm27_improper(lfailure,ffname,nimprotypes)
c     ******************************************************************
c     * Charmm27 improper harmonic torsions                            *
c     *                                                                *
c     * split from ffcharmm27 03-05-2006 by M.G. Martin                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nimprotypes
c     --- local variables
      integer type,isame
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)
      integer nislocal 
      character*(FFNAMELEN) ninlocal 
      dimension ninlocal(MAXSAMEIMP,MAXIMPNAME)
      double precision onepi

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

      call twh_constant_pi(GLB_GET,onepi)

c     --- initialize type count
      type = 0

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (100.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- OBL - X - X - CL [par_all27_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'OBL'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CL'
c     --- NY - CA - CY - CPT [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NY'
      ninlocal(isame,2) = 'CA'
      ninlocal(isame,3) = 'CY'
      ninlocal(isame,4) = 'CPT'
c     --- OB - X - X - CD [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'OB'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CD'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (3.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HEL2 - HEL2 - CEL2 - CEL2 [par_all27_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HEL2'
      ninlocal(isame,2) = 'HEL2'
      ninlocal(isame,3) = 'CEL2'
      ninlocal(isame,4) = 'CEL2'
c     --- HE2 - HE2 - CE2 - CE2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HE2'
      ninlocal(isame,2) = 'HE2'
      ninlocal(isame,3) = 'CE2'
      ninlocal(isame,4) = 'CE2'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (96.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- OCL - X - X - CL [par_all27_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'OCL'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CL'
c     --- OC - X - X - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'OC'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CC'
c     --- CC - X - X - CT1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CC'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CT1'
c     --- CC - X - X - CT2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CC'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CT2'
c     --- CC - X - X - CT3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CC'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CT3'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

c     *** Protein improper groups (if not listed above) 
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (20.80d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CPB - CPA - NPH - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CPB'
      ninlocal(isame,2) = 'CPA'
      ninlocal(isame,3) = 'NPH'
      ninlocal(isame,4) = 'CPA'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (90.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CPB - X - X - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CPB'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'C'
c     --- CT2 - X - X - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CT2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CPB'
c     --- CT3 - X - X - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CT3'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CPB'
c     --- CN1 - X - X - ON1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN1'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'ON1'
c     --- CN1T - X - X - ON1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN1T'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'ON1'
c     --- CN1 - NN2G - CN5G - ON1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN1'
      ninlocal(isame,2) = 'NN2G'
      ninlocal(isame,3) = 'CN5G'
      ninlocal(isame,4) = 'ON1'
c     --- CN1 - NN2U - CN3T - ON1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN1'
      ninlocal(isame,2) = 'NN2U'
      ninlocal(isame,3) = 'CN3T'
      ninlocal(isame,4) = 'ON1'
c     --- CN2 - X - X - NN1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'NN1'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (20.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - C - C - HA [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HA'
      ninlocal(isame,2) = 'C'
      ninlocal(isame,3) = 'C'
      ninlocal(isame,4) = 'HA'
c     --- HA - CPB - C - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HA'
      ninlocal(isame,2) = 'CPB'
      ninlocal(isame,3) = 'C'
      ninlocal(isame,4) = 'C'
c     --- NH1 - X - X - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NH1'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'H'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (20.00d0)*(CON_KCALTOK)
      iclocal(1) = (180.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - HA - C - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HA'
      ninlocal(isame,2) = 'HA'
      ninlocal(isame,3) = 'C'
      ninlocal(isame,4) = 'C'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (29.40d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - CPA - CPA - CPM [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HA'
      ninlocal(isame,2) = 'CPA'
      ninlocal(isame,3) = 'CPA'
      ninlocal(isame,4) = 'CPM'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.50d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HR1 - NR1 - NR2 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HR1'
      ninlocal(isame,2) = 'NR1'
      ninlocal(isame,3) = 'NR2'
      ninlocal(isame,4) = 'CPH2'
c     --- HR1 - NR2 - NR1 - CPH2 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HR1'
      ninlocal(isame,2) = 'NR2'
      ninlocal(isame,3) = 'NR1'
      ninlocal(isame,4) = 'CPH2'
c     --- HR3 - CPH1 - NR1 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HR3'
      ninlocal(isame,2) = 'CPH1'
      ninlocal(isame,3) = 'NR1'
      ninlocal(isame,4) = 'CPH1'
c     --- HR3 - CPH1 - NR2 - CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HR3'
      ninlocal(isame,2) = 'CPH1'
      ninlocal(isame,3) = 'NR2'
      ninlocal(isame,4) = 'CPH1'
c     --- HR3 - NR1 - CPH1- CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HR3'
      ninlocal(isame,2) = 'NR1'
      ninlocal(isame,3) = 'CPH1'
      ninlocal(isame,4) = 'CPH1'
c     --- HR3 - NR2 - CPH1- CPH1 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HR3'
      ninlocal(isame,2) = 'NR2'
      ninlocal(isame,3) = 'CPH1'
      ninlocal(isame,4) = 'CPH1'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- N - C - CP1 - CP3 [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'N'
      ninlocal(isame,2) = 'C'
      ninlocal(isame,3) = 'CP1'
      ninlocal(isame,4) = 'CP3'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (40.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NC2 - X - X - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NC2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'C'
c     --- CN2 - NN3G - NN2G - NN1 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN2'
      ninlocal(isame,2) = 'NN3G'
      ninlocal(isame,3) = 'NN2G'
      ninlocal(isame,4) = 'NN1'
c     --- CN2 - NN3A - CN5 - NN1 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN2'
      ninlocal(isame,2) = 'NN3A'
      ninlocal(isame,3) = 'CN5'
      ninlocal(isame,4) = 'NN1'
c     --- ON1 - X - X - CN1A [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'ON1'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CN1A'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (4.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NH2 - X - X - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NH2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'H'
c     --- HN1 - X - X - NN1 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN1'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'NN1'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (40.60d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NPH - CPA - CPB - CPB [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NPH'
      ninlocal(isame,2) = 'CPA'
      ninlocal(isame,3) = 'CPB'
      ninlocal(isame,4) = 'CPB'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (18.30d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NPH - CPA - CPM - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NPH'
      ninlocal(isame,2) = 'CPA'
      ninlocal(isame,3) = 'CPM'
      ninlocal(isame,4) = 'CPA'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (32.70d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NPH - CPA - CPM - CPA [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NPH'
      ninlocal(isame,2) = 'CPM'
      ninlocal(isame,3) = 'CPB'
      ninlocal(isame,4) = 'CPA'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (1.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HR3 - CPH1 - NR3 - CPH1[par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HR3'
      ninlocal(isame,2) = 'CPH1'
      ninlocal(isame,3) = 'NR3'
      ninlocal(isame,4) = 'CPH1'
c     --- HN2 - X - X - NN2 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'NN2'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.45d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - CPH1 - CPH2 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NR1'
      ninlocal(isame,2) = 'CPH1'
      ninlocal(isame,3) = 'CPH2'
      ninlocal(isame,4) = 'H'
c     --- NR1 - CPH2 - CPH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NR1'
      ninlocal(isame,2) = 'CPH2'
      ninlocal(isame,3) = 'CPH1'
      ninlocal(isame,4) = 'H'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (1.20d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NR3 - CPH1 - CPH2 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NR3'
      ninlocal(isame,2) = 'CPH1'
      ninlocal(isame,3) = 'CPH2'
      ninlocal(isame,4) = 'H'
c     --- NR3 - CPH2 - CPH1 - H [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NR3'
      ninlocal(isame,2) = 'CPH2'
      ninlocal(isame,3) = 'CPH1'
      ninlocal(isame,4) = 'H'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (45.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- O - CP1 - NH2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'CP1'
      ninlocal(isame,3) = 'NH2'
      ninlocal(isame,4) = 'CC'
c     --- O - CT1 - NH2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'CT1'
      ninlocal(isame,3) = 'NH2'
      ninlocal(isame,4) = 'CC'
c     --- O - CT2 - NH2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'CT2'
      ninlocal(isame,3) = 'NH2'
      ninlocal(isame,4) = 'CC'
c     --- O - CT3 - NH2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'CT3'
      ninlocal(isame,3) = 'NH2'
      ninlocal(isame,4) = 'CC'
c     --- O - HA - NH2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'HA'
      ninlocal(isame,3) = 'NH2'
      ninlocal(isame,4) = 'CC'
c     --- O - NH2 - CP1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'NH2'
      ninlocal(isame,3) = 'CP1'
      ninlocal(isame,4) = 'CC'
c     --- O - NH2 - CT1 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'NH2'
      ninlocal(isame,3) = 'CT1'
      ninlocal(isame,4) = 'CC'
c     --- O - NH2 - CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'NH2'
      ninlocal(isame,3) = 'CT2'
      ninlocal(isame,4) = 'CC'
c     --- O - NH2 - CT3 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'NH2'
      ninlocal(isame,3) = 'CT3'
      ninlocal(isame,4) = 'CC'
c     --- O - NH2 - HA - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'NH2'
      ninlocal(isame,3) = 'HA'
      ninlocal(isame,4) = 'CC'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (120.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- O - N - CT2 - CC [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'N'
      ninlocal(isame,3) = 'CT2'
      ninlocal(isame,4) = 'CC'
c     --- O - X - X - C [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'O'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'C'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

c     *** Nucleic acid improper terms (if not listed above) 
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (7.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NN2B - CN4 - CN5 - HN2 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NN2B'
      ninlocal(isame,2) = 'CN4'
      ninlocal(isame,3) = 'CN5'
      ninlocal(isame,4) = 'HN2'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.80d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NN2G - CN4 - CN1 - HN2 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NN2G'
      ninlocal(isame,2) = 'CN4'
      ninlocal(isame,3) = 'CN1'
      ninlocal(isame,4) = 'HN2'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (6.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NN1 - CN2 - HN1 - HN1 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NN1'
      ninlocal(isame,2) = 'CN2'
      ninlocal(isame,3) = 'HN1'
      ninlocal(isame,4) = 'HN1'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (110.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1T - NN2B - NN2U - ON1 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN1T'
      ninlocal(isame,2) = 'NN2B'
      ninlocal(isame,3) = 'NN2U'
      ninlocal(isame,4) = 'ON1'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (80.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - X - X - ON1C [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN1'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'ON1C'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (60.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CN2 - NN3 - CN3 - NN1 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN2'
      ninlocal(isame,2) = 'NN3'
      ninlocal(isame,3) = 'CN3'
      ninlocal(isame,4) = 'NN1'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (39.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CN4 - NN2G - NN3I - HN3 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN4'
      ninlocal(isame,2) = 'NN2G'
      ninlocal(isame,3) = 'NN3I'
      ninlocal(isame,4) = 'HN3'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (14.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CN9 - X - X - CN3T [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN9'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CN3T'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (15.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CN3 - CN3C - CN8 - HN6 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'CN3'
      ninlocal(isame,2) = 'CN3C'
      ninlocal(isame,3) = 'CN8'
      ninlocal(isame,4) = 'HN6'
c     --- HN3B - X - X - CN3 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN3B'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CN3'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (13.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HN3B - X - X - CN3A [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN3B'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CN3A'
c     --- HN3B - X - X - CN3B [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN3B'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CN3B'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (50.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HN2 - CN3 - CN3B - NN2 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN2'
      ninlocal(isame,2) = 'CN3'
      ninlocal(isame,3) = 'CN3B'
      ninlocal(isame,4) = 'NN2'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (-5.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HN1 - HN1 - CN1A - NN1 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN1'
      ninlocal(isame,2) = 'HN1'
      ninlocal(isame,3) = 'CN1A'
      ninlocal(isame,4) = 'NN1'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (53.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HN3 - X - X - CN3C [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN3'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CN3C'
c     --- HN6 - X - X - CN3C [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN6'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'CN3C'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (18.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- HN8 - CN3 - CN3 - CN8 [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'HN8'
      ninlocal(isame,2) = 'CN3'
      ninlocal(isame,3) = 'CN3'
      ninlocal(isame,4) = 'CN8'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.60d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - CPH1 - CPH2 - CN7B [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NR1'
      ninlocal(isame,2) = 'CPH1'
      ninlocal(isame,3) = 'CPH2'
      ninlocal(isame,4) = 'CN7B'
c     --- NR1 - CPH2 - CPH1 - CN7B [par_all27_na.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NR1'
      ninlocal(isame,2) = 'CPH2'
      ninlocal(isame,3) = 'CPH1'
      ninlocal(isame,4) = 'CN7B'
c     --- set nimpsame
      nislocal = isame
        call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (137.40d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NPH - CPA - CPA - FE [par_all27_prot_lipid.inp]
      isame = isame + 1
      ninlocal(isame,1) = 'NPH'
      ninlocal(isame,2) = 'CPA'
      ninlocal(isame,3) = 'CPA'
      ninlocal(isame,4) = 'FE'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return
      
      return
      end
