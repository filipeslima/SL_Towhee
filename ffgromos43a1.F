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
      subroutine twh_ffgromos43a1(lfailure,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the gromos43a1   *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-20-2002 by M.G. Martin  *
c     * last modified 07-21-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_NBCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
     &        ,nimprotypes
c     --- local variables
      logical lrecompute
      character*(FFNAMELEN) ffname
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      integer i,j,ii,jj,ij,type
      integer jlo,jhi,itrial,select,isel,jsel
      dimension select(44,44)
      double precision csixi,csixj,ctwvi,ctwvj,csixt,ctwvt
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)
      integer nislocal 
      double precision onepi
c     --- set up the Gromos selection matrix
      data (select(1,i),i=1,44) /1,1,2,2,2,2,2,2,2,2,10*1,2,2,2
     &     ,2,2,2,2,1*1,1,1,1,1*1,2,2,1,5*1,1,2*1,2/
      data (select(2,i),i=1,44) /1,1,2,2,2,2,3,2,3,3,10*1,3,3,3
     &     ,3,3,3,3,1*1,1,1,1,1*1,2,3,1,5*1,1,2*1,3/
      data (select(3,i),i=1,44) /2,2,2,2,2,2,2,2,2,2,10*1,2,2,2
     &     ,2,2,2,2,1*1,2,2,2,1*1,2,2,2,5*1,2,2*1,2/
      data (select(4,i),i=1,44) /2,2,2,2,2,2,2,2,2,2,10*1,2,2,2
     &     ,2,2,2,2,1*1,2,2,2,1*1,2,2,2,5*1,2,2*1,2/
      data (select(5,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(6,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(7,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(8,i),i=1,44) /2,2,2,2,2,2,2,2,2,2,10*1,2,2,2
     &     ,2,2,2,2,1*1,2,2,2,1*1,2,2,2,5*1,2,2*1,2/
      data (select(9,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(10,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(11,i),i=1,44) /44*1/ 
      data (select(12,i),i=1,44) /44*1/ 
      data (select(13,i),i=1,44) /44*1/ 
      data (select(14,i),i=1,44) /44*1/ 
      data (select(15,i),i=1,44) /44*1/ 
      data (select(16,i),i=1,44) /44*1/ 
      data (select(17,i),i=1,44) /44*1/ 
      data (select(18,i),i=1,44) /44*1/ 
      data (select(19,i),i=1,44) /44*1/ 
      data (select(20,i),i=1,44) /44*1/ 
      data (select(21,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(22,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(23,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(24,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(25,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(26,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(27,i),i=1,44) /44*1/
      data (select(28,i),i=1,44) /44*1/
      data (select(29,i),i=1,44) /1,1,2,2,2,2,2,2,2,2,10*1,2,2,2
     &     ,2,2,2,2,1*1,1,1,1,1*1,2,2,1,5*1,1,2*1,2/
      data (select(30,i),i=1,44) /44*1/
      data (select(31,i),i=1,44) /44*1/
      data (select(32,i),i=1,44) /44*1/
      data (select(33,i),i=1,44) /2,2,2,2,2,2,2,2,2,2,10*1,2,2,2
     &     ,2,2,2,2,1*1,2,2,2,1*1,2,2,2,5*1,2,2*1,2/
      data (select(34,i),i=1,44) /2,2,2,2,1,1,1,2,1,1,10*1,1,1,1
     &     ,1,1,1,1,1*1,2,2,2,1*1,2,1,2,5*1,2,2*1,1/
      data (select(35,i),i=1,44) /1,1,2,2,2,2,3,2,3,3,10*1,3,3,3
     &     ,3,3,3,3,1*1,1,1,1,1*1,2,3,1,5*1,1,2*1,3/
      data (select(36,i),i=1,44) /36*1,0,0,5*1,1/
      data (select(37,i),i=1,44) /35*1,0,1,0,5*1,1/
      data (select(38,i),i=1,44) /35*1,0,0,6*1,1/
      data (select(39,i),i=1,44) /39*1,0,0,2*1,1/
      data (select(40,i),i=1,44) /38*1,0,1,0,2*1,1/
      data (select(41,i),i=1,44) /1,1,2,2,2,2,2,2,2,2,10*1,2,2,2
     &     ,2,2,2,2,1*1,1,1,1,1*1,2,2,1,4*1,0,1,2*1,2/
      data (select(42,i),i=1,44) /44*1/
      data (select(43,i),i=1,44) /44*1/
      data (select(44,i),i=1,44) /44*1/
c     --- end of the Gromos selection matrix

      ffname = 'Gromos43A1'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     Note: The builder will not work with Gromos as it does not have
c     an easily discerned naming convention.  All Naming convention
c     automation has been removed to reflect this.

c     Gromos force field has special mixing rules that are 
c     not easily combined with other force fields.  Thus if you 
c     are using the Gromos mixrule you have to use only the Gromos
c     force field

c     **************************************************************
c     * 1-44 GROMOS force field Gromos 96 43A1                     *
c     **************************************************************
c     the type numbers now are the same as the gromos manual numbers
c     parameters are listed using gromos units (except for Si which 
c     is given number 44 so it can have a different mass from P
c     --- csix: sqrt[ kJ nm^6 / mol ]
c     --- ctwelve: 10^(-3) sqrt[ kJ nm^12 / mol ]
c     selection matrix:
c     --- Gromos has up to three different c12half parameters and 
c     --- this matrix contains the information of which one to use 
c     --- for the different combinations of atom types
c     --- this matrix is set up at the top of this subroutine 
c     references:
c     [Gromos manual] W.F. van Gunsteren, S.R. Billeter, 
c     --------------- A.A. Eising, P.H. Hunenberger, P. Kruger,
c     --------------- A.E. Mark, W.R.P. Scott, I.G. Tironi
c     --------------- `Biomolecular Simulation: The GROMOS96 manual
c     --------------- and user guide.' VdF: Hochschulverlag AG an 
c     --------------- der ETH Zurich and BIOMOS b.v, Zurich, 
c     --------------- Groningen, 1996 ISBN 3 7281 2422 2.
c     --------------- Good luck trying to find this without buying
c     --------------- Gromos as I couln't find anyone who would sell
c     --------------- me a copy.  My copy was kindly provided by 
c     --------------- Gromos developer Lukas D. Schuler (MG Martin)
      
c     (O) carbonyl oxygen C=O [Gromos manual]
      type = 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04756d0
      nbclocal(2) = 0.8611d0
      nbclocal(3) = 1.125d0
      masslocal = 15.999d0
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (OM) carboxyl oxygen C-O [Gromos manual]
      type = 2
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04756d0
      nbclocal(2) = 0.8611d0
      nbclocal(3) = 1.841d0
      nbclocal(4) = 3.068d0
      masslocal = 15.999d0
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OM'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (OA) hydroxyl, sugar or ester Oxygen [Gromos manual]
c     --- broken into 2 groups for bonds (MGM 06-22-2001)
c     --- OAa: aliphatic 
c     --- OAs: sugar
      type = 3
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04756d0
      nbclocal(2) = 1.125d0
      nbclocal(3) = 1.227d0
      masslocal = 15.999d0
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OA'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (OW) water oxygen [Gromos manual]
      type = 4
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.05116d0
      nbclocal(2) = 1.544d0
      nbclocal(3) = 1.623d0
      masslocal = 15.999d0
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OW'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (N) peptide nitrogen NH [Gromos manual]
      type = 5
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04936d0
      nbclocal(2) = 1.301d0
      nbclocal(3) = 1.943d0
      masslocal = 14.007d0
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (NT) terminal nitrogen NH2 [Gromos manual]
      type = 6
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04936d0
      nbclocal(2) = 1.301d0
      nbclocal(3) = 2.250d0
      masslocal = 14.007d0
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='NT'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (NL) terminal nitrogen NH3 [Gromos manual]
      type = 7
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04936d0
      nbclocal(2) = 1.301d0
      nbclocal(3) = 3.068d0
      masslocal = 14.007d0
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='NL'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (NR) aromatic nitrogen [Gromos manual]
      type = 8
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04936d0
      nbclocal(2) = 1.301d0
      nbclocal(3) = 1.841d0
      masslocal = 14.007d0
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='NR'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (NZ) Arg NH (NH2) [Gromos manual]
      type = 9
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04936d0
      nbclocal(2) = 1.301d0
      nbclocal(3) = 2.148d0
      masslocal = 14.007d0
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='NZ'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (NE) Arg NE (NH) [Gromos manual]
      type = 10
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04936d0
      nbclocal(2) = 1.301d0
      nbclocal(3) = 1.984d0
      masslocal = 14.007d0
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='NE'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (C) bare carbon [Gromos manual]
      type = 11
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04838d0
      nbclocal(2) = 1.837d0
      masslocal = 12.011d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- NOTE: For atom types 12 to 16 (except 15) the one-four 
c     --- interactions are different from the normal LJ interactions
c     --- These one-four interactions are stored as follows
c     --- one-four c6_half is in ctwelve(*,2)
c     --- one-four c12_half is in ctwelve(*,3)
      
c     (CH1) aliphatic or sugar CH group [Gromos manual]
      type = 12
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.06148d0
      nbclocal(2) = 3.373d0
c     --- special, store one-four csix and ctwelve in 3 and 4
      nbclocal(3) = 0.05396d0
      nbclocal(4) = 1.933d0
      masslocal = 13.0189d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH1'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (CH2) aliphatic or sugar CH2 group Gromos [Gromos manual]
      type = 13
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.08429d0
      nbclocal(2) = 5.077d0
c     --- special, store one-four csix and ctwelve in 3 and 4
      nbclocal(3) = 0.06873d0
      nbclocal(4) = 2.667d0
      masslocal = 14.0268d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH2'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (CH3) aliphatic CH3 group [Gromos manual]
      type = 14
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.09958d0
      nbclocal(2) = 5.794d0
c     --- special, store one-four csix and ctwelve in 3 and 4
      nbclocal(3) = 0.08278d0
      nbclocal(4) = 3.473d0
      masslocal = 15.0347d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (CH4) methane [Gromos manual]
      type = 15
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.1148d0
      nbclocal(2) = 5.862d0
      masslocal = 16.043d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH4'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (CR1) Aromatic CH group [Gromos manual]
      type = 16
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.07425d0
      nbclocal(2) = 3.888d0
c     --- special, store one-four csix and ctwelve in 3 and 4
      nbclocal(3) = 0.07435d0
      nbclocal(4) = 2.886d0
      masslocal = 13.0189d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CR1'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (HC) hydrogen bound to carbon [Gromos manual]
c     --- used for hydrogens bonded to aromatic carbons as aliphatic
c     --- are lumped into united-atoms
      type = 17
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.00928d0
      nbclocal(2) = 0.123d0
      masslocal = 1.0079d0
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HC'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (H) hydrogen not bound to carbon [Gromos manual]
      type = 18
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      masslocal = 1.0079d0
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (DUM) dummy atom [Gromos manual]
      type = 19
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      elelocal = '--'
      masslocal = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='DUM'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (S) sulphur [Gromos manual]
      type = 20
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.09992d0
      nbclocal(2) = 3.616d0
      masslocal = 32.066d0
      elelocal = ' S'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='S'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Cu1+) copper with charge +1 [Gromos manual]
      type = 21
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.02045d0
      nbclocal(2) = 0.07159d0
      nbclocal(3) = 0.2250d0
      masslocal = 63.546d0
      elelocal = 'Cu'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cu1+'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'ion'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Cu2+) copper with charge +2 [Gromos manual]
      type = 22
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.02045d0
      nbclocal(2) = 0.07159d0
      nbclocal(3) = 0.4091d0
      masslocal = 63.546
      elelocal = 'Cu'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cu2+'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'ion'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Fe) iron in heme [Gromos manual]
      type = 23
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      masslocal = 55.847d0
      elelocal = 'Fe'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Fe'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'heme'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Zn2+) zinc with charge 2+ [Gromos manual]
      type = 24
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.02045d0
      nbclocal(2) = 0.09716d0
      nbclocal(3) = 0.09716d0
      masslocal = 65.39d0
      elelocal = 'Zn'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Zn2+'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'ion'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Mg2+) magnesium with charge 2+ [Gromos manual]
      type = 25
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.018080d0
      nbclocal(2) = 0.05838d0
      nbclocal(3) = 0.05838d0
      masslocal = 24.305d0
      elelocal = 'Mg'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Mg2+'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'ion'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Ca2+) calcium with charge 2+ [Gromos manual]
      type = 26
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.03170d0
      nbclocal(2) = 0.7057d0
      nbclocal(3) = 0.7057d0
      masslocal = 24.305d0
      elelocal = 'Ca'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ca2+'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'ion'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (P) phosphorus  [Gromos manual]
c     --- NOTE: Silicon has the same parameters but a different mass
c     --- it is stored as type 144
      type = 27
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.1214d0
      nbclocal(2) = 4.711d0
      masslocal = 30.974d0
      elelocal = ' P'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='P'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Ar) argon [Gromos manual]
      type = 28
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.07915d0
      nbclocal(2) = 3.138d0
      masslocal = 39.948d0
      elelocal = 'Ar'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ar'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (F) fluorine (non-ionic) [Gromos manual]
      type = 29
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.03432d0
      nbclocal(2) = 0.8722d0
      nbclocal(3) = 1.227d0
      masslocal = 18.998d0
      elelocal = ' F'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='F'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Cl) chlorine (non-ionic)  [Gromos manual]
      type = 30
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.09362d0
      nbclocal(2) = 3.911d0
      masslocal = 18.998d0
      elelocal = 'Cl'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cl'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Br) bromine (non-ionic) [Gromos manual]
      type = 31
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.03434d0
      nbclocal(2) = 8.092d0
      masslocal = 35.453d0
      elelocal = 'Br'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Br'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (CMet) Ch3 group in methanol [Gromos manual]
      type = 32
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.09421d0
      nbclocal(2) = 4.5665d0
      masslocal = 15.0347d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CMet'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (OMet) Oxygen in methanol [Gromos manual]
      type = 33
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.04756d0
      nbclocal(2) = 1.125d0
      nbclocal(3) = 1.227d0
      masslocal = 15.999d0
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OMet'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Na+) Sodium with charge +1 [Gromos manual]
      type = 34
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.008489d0
      nbclocal(2) = 0.1450d0
      nbclocal(3) = 0.1450d0
      masslocal = 22.990d0
      elelocal = 'Na'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Na+'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'ion'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Cl-) Chlorine with charge -1 [Gromos manual]
      type = 35
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.1175d0
      nbclocal(2) = 10.340d0
      nbclocal(3) = 10.340d0
      nbclocal(4) = 10.340d0
      masslocal = 35.453d0
      elelocal = 'Cl'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cl-'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'ion'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (CChl) Carbon in chloroform [Gromos manual]
      type = 36
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.051292d0
      nbclocal(2) = 2.0160d0
      masslocal = 12.011d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CChl'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (ClChl) Chlorine in chloroform [Gromos manual]
      type = 37
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.091141d0
      nbclocal(2) = 3.7101d0
      masslocal = 35.453d0
      elelocal = 'Cl'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='ClChl'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (HChl) hydrogen in chloroform [Gromos manual]
      type = 38
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.00614d0
      nbclocal(2) = 0.065574d0
      masslocal = 1.0079d0
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HChl'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (SDmso) Sulphur in DMSO [Gromos manual]
      type = 39
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.10277d0
      nbclocal(2) = 4.6366d0
      masslocal = 32.06d0
      elelocal = ' S'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='SDmso'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (CDmso) CH3 group in DMSO [Gromos manual]
      type = 40
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.095139d0
      nbclocal(2) = 4.6645d0
      masslocal = 15.0347d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CDmso'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (ODmso) Oxygen in DMSO [Gromos manual]
      type = 41
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.047652d0
      nbclocal(2) = 0.86686d0
      nbclocal(3) = 1.125d0
      masslocal = 15.999d0
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='ODmso'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (CCl4) Carbon in carbontetrachloride [Gromos manual]
      type = 42
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.051292d0
      nbclocal(2) = 2.7568d0
      masslocal = 12.011d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CCl4'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (ClCl4) Chlorine in carbontetrachloride [Gromos manual]
      type = 43
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.087201d0
      nbclocal(2) = 3.5732d0
      masslocal = 12.011d0
      elelocal = 'Cl'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='ClCl4'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (Si) silicon [Gromos manual]
c     --- this has the same parameters as P(type127) but a different
c     --- element and therefore a different mass
      type = 44
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.1214d0
      nbclocal(2) = 4.711d0
      masslocal = 28.086d0
      elelocal = 'Si'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Si'
      call twh_nbname(GLB_SET,type,nbnlocal)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     *****************************************************************
c     * use the general nonbond terms plus the gromos selection       *
c     * matrix to populate all of the cross terms as they need to be  *
c     * explicitly listed for Grommos                                 *
c     *****************************************************************

c     --- initialize
      csixi = 0.0d0
      csixj = 0.0d0
      csixt = 0.0d0
      ctwvi = 0.0d0
      ctwvj = 0.0d0
      ctwvt = 0.0d0

      do itrial = 1,2
         do i = 1, ntypes
            if ( itrial .eq. 1 ) then
c               --- we are setting up all of the dissimilar cross terms
               jlo = i+1
               jhi = ntypes
            else
c              --- we are setting up all of the self terms
               jlo = i
               jhi = i
            endif
            do j = jlo,jhi
               ii = twh_getnbtype(i,i)
               jj = twh_getnbtype(j,j)
               ij = twh_getnbtype(i,j)
c              --- which currently range only up to type 44
               if ( i .le. 44 .and. j .le. 44 ) then
c                 --- compute the gromos selection matrix
                  isel = select(i,j)
                  jsel = select(j,i)
                  if ( isel .eq. 0 .or. jsel .eq. 0 ) then
c                    --- special cross term: consult if list
                     if ( (i .eq. 36 .and. j .eq. 37) .or.
     &                    (i .eq. 37 .and. j .eq. 36) ) then
                        csixt = 4.6754d-3
                        ctwvt = 7.4813d0
                     elseif ( (i .eq. 36 .and. j .eq. 38) .or.
     &                       (i .eq. 38 .and. j .eq. 36) ) then
                        csixt = 0.3622d-3
                        ctwvt = 0.1745d0
                     elseif ( (i .eq. 37 .and. j .eq. 38) .or.
     &                       (i .eq. 38 .and. j .eq. 37)) then
                        csixt = 0.6493d-3
                        ctwvt = 0.3266d0
                     elseif ( (i .eq. 39 .and. j .eq. 40) .or.
     &                       (i .eq. 40 .and. j .eq. 39)) then
                        csixt = 9.7827d-3
                        ctwvt = 21.6523d0
                     elseif ( (i .eq. 39 .and. j .eq. 41) .or.
     &                       (i .eq. 41 .and. j .eq. 39) ) then
                        csixt = 5.2442d-3
                        ctwvt = 4.6094d0
                     elseif ( (i .eq. 40 .and. j .eq. 41) .or.
     &                       (i .eq. 41 .and. j .eq. 40) ) then
                        csixt = 4.9187d-3
                        ctwvt = 4.7597d0
                     endif
                  else
c                    --- use gromos selection matrix to determine the
c                    --- parameters
                     csixi = twh_get_nbcoeff(1,ii)
                     csixj = twh_get_nbcoeff(1,jj)
                     csixt = csixi*csixj
                     ctwvi = twh_get_nbcoeff((1+isel),ii)
                     ctwvj = twh_get_nbcoeff((1+jsel),jj)
                     ctwvt = ctwvi*ctwvj
                  endif
                     
                  if ( ctwvt .eq. 0.0d0 ) then
                     nbclocal(1) = 0.0d0
                     nbclocal(2) = 0.0d0
                     nbclocal(3) = 0.0d0
                     nbclocal(4) = 0.0d0
                  else
c                    --- convert to standard K and A units
                     csixt = csixt*(CON_KJTOK)*(1d6)
                     ctwvt = ctwvt*(1d-6)*(CON_KJTOK)*(1d12)
c                    --- convert to sigma and epsilon
                     nbclocal(1) = (ctwvt/csixt)**(1.0d0/6.0d0)
                     nbclocal(2) = 0.25d0*(csixt*csixt)/ctwvt

                     lrecompute = .false.
c                    --- specify the one-four interactions
                     if ( i .eq. 12 .or. i .eq. 13  .or. i .eq. 14
     &                    .or. i .eq. 16) then
c                       --- special one-four interaction parameters
                        lrecompute = .true.
                        csixi = twh_get_nbcoeff(3,ii)
                        ctwvi = twh_get_nbcoeff(4,ii)
                     endif
                     if ( j .eq. 12 .or. j .eq. 13  .or. j .eq. 14
     &                    .or. j .eq. 16) then
c                       --- special one-four interaction parameters
                        lrecompute = .true.
                        csixj = twh_get_nbcoeff(3,jj)
                        ctwvj = twh_get_nbcoeff(4,jj)
                     endif
                     if ( lrecompute ) then
                        csixt = csixi*csixj*(CON_KJTOK)*(1d6)
                        ctwvt = ctwvi*ctwvj*(1d-6)*(CON_KJTOK)*(1d12)
                     endif
c                    --- convert to sigma and epsilon
                     nbclocal(3) = (ctwvt/csixt)**(1.0d0/6.0d0)
                     nbclocal(4) = 0.25d0*(csixt*csixt)/ctwvt
                  endif
                  call twh_set_nbcoeff(ij,nbclocal)
               endif
            enddo
         enddo
      enddo

c     ******************************************************************
c     * 1-47 GROMOS-43A1 bond potentials                               *
c     ******************************************************************
c     --- Gromos Manual, Gromos96 43A1 force field
c     --- the type are the same as the Gromos bond types
c     --- This is from Table 2.5.2.1 in the gromos manual

c     --- H - OA
      type = 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.0d0
      vclocal(1) = (15.7d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- H - N (all)
      type = 2
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.0d0
      vclocal(1) = (18.7d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- HC - C
      type = 3
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.09d0
      vclocal(1) = (12.3d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C - O
      type = 4
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.23d0
      vclocal(1) = (16.6d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C - OM
      type = 5
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.25d0
      vclocal(1) = (13.4d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CR1 - NR (6-ring)
      type = 6
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.32d0
      vclocal(1) = (12.0d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- H - S
      type = 7
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.33d0
      vclocal(1) = (8.87d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C - NT,NL
      type = 8
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.33d0
      vclocal(1) = (10.6d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C,CR1 - N,NR,CR1,C (peptide, 5-ring)
      type = 9
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.33d0
      vclocal(1) = (11.8d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C - N,NZ,NE
      type = 10
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.34d0
      vclocal(1) = (10.5d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C - NR (no H) (6-ring)
      type = 11
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.34d0
      vclocal(1) = (11.7d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C - OA (both kinds of OA ?)
      type = 12
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.36d0
      vclocal(1) = (10.2d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C - NR (heme)
      type = 13
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.38d0
      vclocal(1) = (11.0d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CH2 - C,CR1 (6-ring)
      type = 14
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.39d0
      vclocal(1) = (8.66d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C,CR1 - CH2,C,CR1 (6-ring)
      type = 15
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.39d0
      vclocal(1) = (10.8d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C,CR1,CH2 - NR (6-ring)
      type = 16
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.40d0
      vclocal(1) = (8.54d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CHn - OA (non-sugar)
      type = 17
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.43d0
      vclocal(1) = (8.18d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CHn - OM
      type = 18 
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.43d0
      vclocal(1) = (9.21d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CHn - OA (sugar)
      type = 19
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.435d0
      vclocal(1) = (6.10d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CHn - N,NT,NL,NZ,NE
      type = 20
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.47d0
      vclocal(1) = (8.71d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CHn - NR (5-ring)
      type = 21
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.48d0
      vclocal(1) = (5.73d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CHn - NR (6-ring)
      type = 22
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.48d0
      vclocal(1) = (7.64d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- O,OM - P
      type = 23
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.48d0
      vclocal(1) = (8.60d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- O - S
      type = 24
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.50d0
      vclocal(1) = (8.37d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CHn - CHn (sugar)
      type = 25
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.52d0
      vclocal(1) = (5.43d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- C,CHn - C,CHn (non-sugar)
      type = 26
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.53d0
      vclocal(1) = (7.15d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- OA - P
      type = 27
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.61d0
      vclocal(1) = (4.84d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- OA - Si
      type = 28
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.63d0
      vclocal(1) = (4.72d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CH3 - S
      type = 29
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.78d0
      vclocal(1) = (5.94d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CH2 - S
      type = 30
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.83d0
      vclocal(1) = (5.62d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CH1 - SI
      type = 31
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.87d0
      vclocal(1) = (3.59d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- NR - FE
      type = 32
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.98d0
      vclocal(1) = (0.640d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- S - S
      type = 33
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 2.04d0
      vclocal(1) = (5.03d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- NR (heme) - FE
      type = 34
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 2.00d0
      vclocal(1) = (0.628d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- HW - OW
      type = 35
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.00d0
      vclocal(1) = (23.2d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- HChl - CChl
      type = 36
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.10d0
      vclocal(1) = (12.1d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CChl - ClChl
      type = 37
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.758d0
      vclocal(1) = (8.12d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- ODmso - SDmso
      type = 38
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.53d0
      vclocal(1) = (8.04d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- SDmso - CDmso
      type = 39
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.95d0
      vclocal(1) = (4.95d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CCl4 - ClCl4
      type = 40
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.76d0
      vclocal(1) = (8.10d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- HW - HW 
      type = 41
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.63299d0
      vclocal(1) = (8.71d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- HChl - CLChl
      type = 42
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 2.33839d0
      vclocal(1) = (2.68d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- ClChl - ClChl
      type = 43
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 2.90283d0
      vclocal(1) = (2.98d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- ODmso - CDmso
      type = 44
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 2.80412d0
      vclocal(1) = (2.39d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- CDmso - CDmso
      type = 45
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 2.92993d0
      vclocal(1) = (2.19d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- HMet - CMet (I have no idea what this means)
      type = 46
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 1.98842d0
      vclocal(1) = (3.97d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

c     --- ClCl4 - ClCl4
      type = 47
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_GROMOS)
      vclocal(0) = 2.87407d0
      vclocal(1) = (3.04d6)*(0.25d0)*(CON_KJTOK)*(1.0d-4)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * 01-46 GROMOS ANGLES harmonic cos                               *
c     ******************************************************************
c     from the Gromos manual
c     same numbering system as in Gromos

c     --- NR(heme) - FE - NR(heme)
      type = 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 90.0d0
      bclocal(1) = (420.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H - S - CH2
      type = 2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 96.0d0
      bclocal(1) = (405.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CH2 - S - CH3
      type = 3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 100.0d0
      bclocal(1) = (475.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- OA - P - OA
      type = 4
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 103.0d0
      bclocal(1) = (420.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CH2 - S - S
      type = 5
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 104.0d0
      bclocal(1) = (490.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- NR - C - CR1 (5-ring)
      type = 6
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 108.0d0
      bclocal(1) = (465.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CHn - CHn - CHn,NR (6-ring) (sugar)
      type = 7
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 109.5d0
      bclocal(1) = (285.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      
c     --- CHn,OA - CHn - OA,NR(5-ring) (sugar)
      type = 8
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 109.5d0
      bclocal(1) = (320.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H - NL,NT - H : CHn - OA - CHn (sugar)
      type = 9
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 109.5d0
      bclocal(1) = (380.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H - NL - C,CHn : H - NT - CHn
      type = 10
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 109.5d0
      bclocal(1) = (425.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- X - OA,SI - X
      type = 11
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 109.5d0
      bclocal(1) = (450.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CHn,C - CHn - C,CHn,OA,OM,N,NE (note similar bond pattern 264)
      type = 12
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 109.5d0
      bclocal(1) = (520.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- OM - P - OA
      type = 13
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 109.6d0
      bclocal(1) = (450.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CHn - CHn - C,CHn,OA,NR,NT,NL (note similar bond pattern 262)
      type = 14
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 111.0d0
      bclocal(1) = (530.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CHn - CH2 - S
      type = 15
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 113.0d0
      bclocal(1) = (545.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- NR(heme) - FE - NR
      type = 16
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 115.0d0
      bclocal(1) = (50.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H - N - CHn
      type = 17
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 115.0d0
      bclocal(1) = (460.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CHn,C - C -OA,N,NT,NL
      type = 18
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 115.0d0
      bclocal(1) = (610.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H -NE - CH2
      type = 19
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 116.0d0
      bclocal(1) = (465.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CH2 - N - CH1
      type = 20
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 116.0d0
      bclocal(1) = (620.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CH3 - N - C : CHn - C - OM
      type = 21
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 117.0d0
      bclocal(1) = (635.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H - NT,NZ,NE - C
      type = 22
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 120.0d0
      bclocal(1) = (390.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H -NT,NZ - H
      type = 23
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 120.0d0
      bclocal(1) = (445.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H - N - CH3,H,HC (6-ring) : H - NT - CHn
      type = 24
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 120.0d0
      bclocal(1) = (505.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- P,SI - OA - CHn,P
      type = 25
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 120.0d0
      bclocal(1) = (530.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- N,C,CR1 (6-ring, no H)
      type = 26
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 120.0d0
      bclocal(1) = (560.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- NZ - C - NZ,NE
      type = 27 
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 120.0d0
      bclocal(1) = (670.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- OM - P - OM
      type = 28
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 120.0d0
      bclocal(1) = (780.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- O - C - CHn,C : CH3 - N 0- CHn
      type = 29
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 121.0d0
      bclocal(1) = (685.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CH1,CH2 - N - C
      type = 30
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 122.0d0
      bclocal(1) = (700.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H - N - C
      type = 31
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 123.0d0
      bclocal(1) = (415.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- O - C - OA,N,NT,NL : C - NE - CH2
      type = 32
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 124.0d0
      bclocal(1) = (730.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- FE - NR - CR1 (5-ring)
      type = 33
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 125.0d0
      bclocal(1) = (375.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- ??
      type = 34
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 125.0d0
      bclocal(1) = (750.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- H,HC - 5-ring
      type = 35
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 126.0d0
      bclocal(1) = (575.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- X(noH) - 5-ring
      type = 36
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 126.0d0
      bclocal(1) = (640.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- OM -C - OM
      type = 37
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 126.0d0
      bclocal(1) = (770.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- 5,6 ring connection
      type = 38
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 132.0d0
      bclocal(1) = (760.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- SI - OA - SI
      type = 39
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 155.0d0
      bclocal(1) = (2215.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- HWat - OWat - HWat
      type = 40
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 109.5d0
      bclocal(1) = (434.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- HChl - CChl - CLChl
      type = 41
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 107.57d0
      bclocal(1) = (484.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CLChl - CChl - CLChl
      type = 42
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 111.30d0
      bclocal(1) = (632.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CDmso - SDmso - CDmso
      type = 43
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 97.4d0
      bclocal(1) = (469.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CDmso - SDmso - ODmso
      type = 44
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 106.75d0
      bclocal(1) = (503.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- HMet - OMet - CMet
      type = 45
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 108.53d0
      bclocal(1) = (443.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

c     --- CLCl4 - CCl4 - CLCl4
      type = 46
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 109.5d0
      bclocal(1) = (618.0d0)*(0.5d0)*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * GROMOS, Charmm, and Amber parameter decoder                    *
c     * Also used for some TraPPE torsions                             *
c     * Sum{1,ntorloop} torcoeff(3(n-1)+1)                             *
c     * *[ 1 + Cos( torcoeff(3(n-1)+2)* theta - torcoeff(3(n-1)+3)mm) ]*
c     *      +  one-four vdw nonbond term                              *
c     *      + torcoeff(-1)*(one-four charge nonbonded term)           *
c     *      this one-four term is only computed if loftor is true     *
c     ******************************************************************

c     ******************************************************************
c     * 1-21 GROMOS standard dihedral                                  *
c     ******************************************************************
c     Gromos 96 43A1 torsional parameters from Gromos manual
c     also includes a 1-4 interaction term (see ffnonbond.F)
c     numbers are the same as the gromos numbers
c     Note: torsion 6 has no contribution from the cosine series, but 
c     is needed in order to get the 1-4 nonbonded terms
c     Note: gromos lists a torsional potential with a +1 or -1 prefactor
c     in front of the cosine term.  This is the same as setting 
c     torcoeff(+3) to 0.0 (for the +1 case) or 180.0 (for the -1 case).

      call twh_constant_pi(GLB_GET,onepi)

c     --- X - C - C - X
      type = 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.86d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p is -1
      tclocal(3) = 180.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - C - OA - X (at ring)
      type = 2
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (7.11d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p is -1
      tclocal(3) = 180.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - C - OA - X (carboxyl)
      type = 3
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (16.7d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p is -1
      tclocal(3) = 180.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - C - N,NT,NE,NZ,NR - X 
      type = 4
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (33.5d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p is -1
      tclocal(3) = 180.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - C - CR1 - X (6-ring)
      type = 5
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (41.8d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p is -1
      tclocal(3) = 180.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - CH1 (sugar) - NR (base) - X
c     --- Null torsion, used for all "extra" torsions in gromos scheme
      type = 6
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,0)
      tclocal(1) = (0.0d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p does not matter (null torsion)
      tclocal(3) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- O - CH1 - CHn - (X != O)
      type = 7
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.418d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- O - CH1 - CHn - O
      type = 8
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.09d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - OA - P - X
      type = 9
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.14d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - S - S - X
      type = 10
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (16.7d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - OA - P - X
      type = 11
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.05d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - CHn - OA (not sugar) - X
      type = 12
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.26d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - CH2 - S - X
      type = 13
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.93d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - C,CHn,SI - NT,NL,OA(sugar) - X
      type = 14
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.77d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- HC - C - S - X
      type = 15
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.18d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- HC - C - C - X
      type = 16
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.44d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - CHn,SI - CHn - X
      type = 17
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.86d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - NR - FE - X
      type = 18
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KJTOK)
      tclocal(2) = 4.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - CHn - N,NE - X
      type = 19
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.0d0)*(CON_KJTOK)
      tclocal(2) = 6.0d0
c     --- gromos p is -1
      tclocal(3) = 180.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - CHn - C,NR(ring),CR1 - X
      type = 20
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.0d0)*(CON_KJTOK)
      tclocal(2) = 6.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

c     --- X - CHn - NT - X
      type = 21
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.77d0)*(CON_KJTOK)
      tclocal(2) = 6.0d0
c     --- gromos p is +1
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * 1-3 Gromos oop Improper Harmonic Torsions                      *
c     * impcoeff(0) * ( phi - impcoeff(1) )^2                          *
c     ******************************************************************
c     Gromos impropers taken from the Gromos manual
c     To convert from the Gromos improper numbering system to this one
c     just add 21 to the Gromos number.
c     As of 8-9-2000 there is no sign checked on the equilibrium angle
c     Thus there is no preference for +35 versus -35 in type 23

c     --- torsional parameter for Gromos planar groups
      type = 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.0510d0)*(CON_KJTOK)
     &     *(0.5d0)*((180.0d0/onepi)**2)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 0
      call twh_nimpsame(GLB_SET,type,nislocal)

c     --- torsional parameter for Gromos tetrahedral centers
      type = 2
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.102d0)*(CON_KJTOK)*(0.5d0)
     &     *((180.0d0/onepi)**2)
      iclocal(1) = (35.26439d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 0
      call twh_nimpsame(GLB_SET,type,nislocal)

c     --- torsional parameter for Gromos heme iron
      type = 3
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.204d0)*(CON_KJTOK)*(0.5d0)
     &     *((180.0d0/onepi)**2)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 0
      call twh_nimpsame(GLB_SET,type,nislocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return

      return
      end
