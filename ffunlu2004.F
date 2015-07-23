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
      subroutine twh_ffunlu2004(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * multiwell force field parameters for SPEAD as described in a   *
c     * series of papers beginning with Unlu et al 2004.               *
c     * The details of this force field were kindly provided by J.R.   *
c     * Elliott Jr. in an email dated 11-25-2005 to M.G. Martin        *
c     *                                                                *
c     * originally written 11-30-2005 by M.G. Martin                   *
c     * last modified 07-21-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
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
      integer isame,type
      double precision center,bondlengtha,bondlengthb
      double precision lawcosa,lawcosb,angcos
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision onepi

      ffname = 'Unlu2004'
c     --- multi-well potential
      call twh_classical_potential(GLB_SET,'Hard 2580 Multistep')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

c     --- set the nonbonded parameters for this force field
      call ffunlu2004_nonbond(type,ffname)

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * square well bonds                                              *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_SQUARE_WELL)
      vclocal(1) = 1.36d0
      vclocal(2) = 1.72d0
      vclocal(3) = 0.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- alkanes [personal communication J.R. Elliott 4-18-2005]
      isame = isame + 1
      vnlocal(isame,1) = 'alkane'
      vnlocal(isame,2) = 'alkane'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_SQUARE_WELL)
      vclocal(1) = 1.29d0
      vclocal(2) = 1.65d0
      vclocal(3) = 0.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- benzene [personal communication J.R. Elliott 4-18-2005]
      isame = isame + 1
      vnlocal(isame,1) = 'benzene'
      vnlocal(isame,2) = 'benzene'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * square well angles                                             *
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)
       
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_SQUARE_WELL)
c     --- set the walls
      bclocal(1) = 2.3356d0
      bclocal(2) = 2.7156d0
      bclocal(3) = 0.0d0
c     --- set a rough equilibrium angle for CBMC
      bondlengtha = 1.54d0
      bondlengthb = 1.54d0
      lawcosa = bondlengtha*bondlengtha + bondlengthb*bondlengthb
      lawcosb = 2.0d0*bondlengtha*bondlengthb
      center = 0.5d0*(bclocal(1)+bclocal(2))
      angcos = (center**2 - lawcosa)/(-lawcosb)
      bclocal(0) = twh_arccos(angcos)*180.0d0/onepi
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- alkane [personal communication J.R. Elliott 4-18-2005]
      isame = isame + 1
      bnlocal(isame,1) = 'alkane'
      bnlocal(isame,2) = 'alkane'
      bnlocal(isame,3) = 'alkane'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_SQUARE_WELL)
c     --- set the walls
      bclocal(1) = 2.36d0
      bclocal(2) = 2.74d0
      bclocal(3) = 0.0d0
c     --- set a rough equilibrium angle for CBMC
      bondlengtha = 1.47d0
      bondlengthb = 1.47d0
      lawcosa = bondlengtha*bondlengtha + bondlengthb*bondlengthb
      lawcosb = 2.0d0*bondlengtha*bondlengthb
      center = 0.5d0*(bclocal(1)+bclocal(2))
      angcos = (center**2 - lawcosa)/(-lawcosb)
      bclocal(0) = twh_arccos(angcos)*180.0d0/onepi
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- benzene [personal communication J.R. Elliott 4-18-2005]
      isame = isame + 1
      bnlocal(isame,1) = 'benzene'
      bnlocal(isame,2) = 'benzene'
      bnlocal(isame,3) = 'benzene'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * square well dihedrals                                          *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,18)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(-1) = 0.0d0
      tclocal(1) = 3.57d0 - 0.15d0
      tclocal(2) = 3.57d0 + 0.15d0
      tclocal(3) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     * - alkane - alkane - * [personal communication J.R. Elliott 
c     4-18-2005]
      tnlocal(1,1) = 'alkane'
      tnlocal(1,2) = 'alkane'
      tnlocal(1,3) = 'alkane'
      tnlocal(1,4) = 'alkane'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,18)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(-1) = 0.0d0
      tclocal(1) = 3.14
      tclocal(2) = 3.425
      tclocal(3) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     * - benzene - benzene - * [personal communication J.R. Elliott 
c     4-18-2005]
      tnlocal(1,1) = 'benzene'
      tnlocal(1,2) = 'benzene'
      tnlocal(1,3) = 'benzene'
      tnlocal(1,4) = 'benzene'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end

      subroutine ffunlu2004_nonbond(type,ffname)
c     ******************************************************************
c     * this subroutine sets up the nonbonded parameters for the unlu  *
c     * et al 2004 force field.                                        *
c     * These parameters are generated via a conversion routine from   *
c     * files provided by J. Richard Elliott                           *
c     *                                                                *
c     * originally written 12-01-2005 by M.G. Martin                   *
c     * last modified 04-17-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer type
      character*(FFNAMELEN) ffname
c     --- local variables
      integer ii
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

c     --- initialize type
      type = 0

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Tang '
      cbnlocal = '    0    0'
      canlocal = '    0    0'
      ctnlocal = '    0    0'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.1000000E+02
      nbclocal(2) =  -0.1560000E+03
      nbclocal(3) =  -0.1040000E+03
      nbclocal(4) =  -0.5200000E+02
      nbclocal(5) =   0.0000000E+00
      nbclocal(6) =   0.1000000E+01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   16.04280
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH4  '
      cbnlocal = '    0    1'
      canlocal = '    0    1'
      ctnlocal = '    0    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3674000E+01
      nbclocal(2) =  -0.1520000E+03
      nbclocal(3) =  -0.1015000E+03
      nbclocal(4) =  -0.5110000E+02
      nbclocal(5) =  -0.6000000E+00
      nbclocal(6) =   0.0000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   39.94800
      elelocal =  'Ar'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ar   '
      cbnlocal = '    0    2'
      canlocal = '    0    2'
      ctnlocal = '    0    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3370000E+01
      nbclocal(2) =  -0.1185000E+03
      nbclocal(3) =  -0.7940000E+02
      nbclocal(4) =  -0.4040000E+02
      nbclocal(5) =  -0.1300000E+01
      nbclocal(6) =   0.0000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.03060
      elelocal =  ' N'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NH3  '
      cbnlocal = '    0    3'
      canlocal = '    0    3'
      ctnlocal = '    0    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3350000E+01
      nbclocal(2) =  -0.3130000E+03
      nbclocal(3) =  -0.2092000E+03
      nbclocal(4) =  -0.1053000E+03
      nbclocal(5) =  -0.1500000E+01
      nbclocal(6) =   0.0000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   34.07580
      elelocal =  ' S'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2S  '
      cbnlocal = '    0    4'
      canlocal = '    0    4'
      ctnlocal = '    0    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3740000E+01
      nbclocal(2) =  -0.2850000E+03
      nbclocal(3) =  -0.1773000E+03
      nbclocal(4) =  -0.9110000E+02
      nbclocal(5) =  -0.7400000E+01
      nbclocal(6) =   0.0000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ' O'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '    0    5'
      canlocal = '    0    5'
      ctnlocal = '    0    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.3500000E+02
      nbclocal(3) =  -0.2850000E+02
      nbclocal(4) =  -0.2210000E+02
      nbclocal(5) =  -0.1560000E+02
      nbclocal(6) =   0.6040000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =    1.00790
      elelocal =  ' H'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H=   '
      cbnlocal = '    0    6'
      canlocal = '    0    6'
      ctnlocal = '    0    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2200000E+01
      nbclocal(2) =   0.2220000E+03
      nbclocal(3) =   0.1478000E+03
      nbclocal(4) =   0.7360000E+02
      nbclocal(5) =  -0.6000000E+00
      nbclocal(6) =   0.3710000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00530
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CO   '
      cbnlocal = '    0    7'
      canlocal = '    0    7'
      ctnlocal = '    0    7'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2830000E+01
      nbclocal(2) =  -0.1520000E+03
      nbclocal(3) =  -0.1015000E+03
      nbclocal(4) =  -0.5110000E+02
      nbclocal(5) =  -0.6000000E+00
      nbclocal(6) =   0.5640000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   38.07180
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'c(=O)'
      cbnlocal = '    0    8'
      canlocal = '    0    8'
      ctnlocal = '    0    8'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2820000E+01
      nbclocal(2) =  -0.1675000E+03
      nbclocal(3) =  -0.1139000E+03
      nbclocal(4) =  -0.6030000E+02
      nbclocal(5) =  -0.6700000E+01
      nbclocal(6) =   0.1160000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   18.99840
      elelocal =  ' F'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F=?  '
      cbnlocal = '    0    9'
      canlocal = '    0    9'
      ctnlocal = '    0    9'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2640000E+01
      nbclocal(2) =  -0.2980000E+02
      nbclocal(3) =  -0.2690000E+02
      nbclocal(4) =  -0.2410000E+02
      nbclocal(5) =  -0.2120000E+02
      nbclocal(6) =   0.7090000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   35.45270
      elelocal =  'Cl'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cl=? '
      cbnlocal = '    0   10'
      canlocal = '    0   10'
      ctnlocal = '    0   10'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3256000E+01
      nbclocal(2) =  -0.1281000E+03
      nbclocal(3) =  -0.9900000E+02
      nbclocal(4) =  -0.6980000E+02
      nbclocal(5) =  -0.4070000E+02
      nbclocal(6) =   0.9940000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =    1.00000
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '2BD- '
      cbnlocal = '    0   11'
      canlocal = '    0   11'
      ctnlocal = '    0   11'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+00
      nbclocal(2) =  -0.2500000E+04
      nbclocal(3) =  -0.2500000E+04
      nbclocal(4) =  -0.2500000E+04
      nbclocal(5) =  -0.2500000E+04
      nbclocal(6) =   0.2700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =    1.00000
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '4BD- '
      cbnlocal = '    0   12'
      canlocal = '    0   12'
      ctnlocal = '    0   12'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+00
      nbclocal(2) =  -0.2500000E+04
      nbclocal(3) =  -0.2500000E+04
      nbclocal(4) =  -0.2500000E+04
      nbclocal(5) =  -0.2500000E+04
      nbclocal(6) =   0.2700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =    1.00000
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '2BA- '
      cbnlocal = '    0   13'
      canlocal = '    0   13'
      ctnlocal = '    0   13'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+00
      nbclocal(2) =  -0.2500000E+04
      nbclocal(3) =  -0.2500000E+04
      nbclocal(4) =  -0.2500000E+04
      nbclocal(5) =  -0.2500000E+04
      nbclocal(6) =   0.2700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =    1.00000
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '4BA- '
      cbnlocal = '    0   14'
      canlocal = '    0   14'
      ctnlocal = '    0   14'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+00
      nbclocal(2) =  -0.2500000E+04
      nbclocal(3) =  -0.2500000E+04
      nbclocal(4) =  -0.2500000E+04
      nbclocal(5) =  -0.2500000E+04
      nbclocal(6) =   0.2700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =    1.00000
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'D-   '
      cbnlocal = '    0   15'
      canlocal = '    0   15'
      ctnlocal = '    0   15'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+00
      nbclocal(2) =  -0.2500000E+04
      nbclocal(3) =  -0.2500000E+04
      nbclocal(4) =  -0.2500000E+04
      nbclocal(5) =  -0.2500000E+04
      nbclocal(6) =   0.2700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3a-'
      cbnlocal = '    1    1'
      canlocal = '    1    1'
      ctnlocal = '    1    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.9090000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4110000E+02
      nbclocal(5) =  -0.1620000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3b-'
      cbnlocal = '    1    2'
      canlocal = '    1    2'
      ctnlocal = '    1    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.5010000E+02
      nbclocal(3) =  -0.4270000E+02
      nbclocal(4) =  -0.3540000E+02
      nbclocal(5) =  -0.2800000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3c-'
      cbnlocal = '    1    3'
      canlocal = '    1    3'
      ctnlocal = '    1    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.3790000E+02
      nbclocal(3) =  -0.3760000E+02
      nbclocal(4) =  -0.3720000E+02
      nbclocal(5) =  -0.3690000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3d-'
      cbnlocal = '    1    4'
      canlocal = '    1    4'
      ctnlocal = '    1    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.4300000E+02
      nbclocal(3) =  -0.4230000E+02
      nbclocal(4) =  -0.4170000E+02
      nbclocal(5) =  -0.4100000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3e-'
      cbnlocal = '    1    5'
      canlocal = '    1    5'
      ctnlocal = '    1    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.6950000E+02
      nbclocal(3) =  -0.5500000E+02
      nbclocal(4) =  -0.4050000E+02
      nbclocal(5) =  -0.2600000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3f-'
      cbnlocal = '    1    6'
      canlocal = '    1    6'
      ctnlocal = '    1    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.1080000E+03
      nbclocal(3) =  -0.7570000E+02
      nbclocal(4) =  -0.4330000E+02
      nbclocal(5) =  -0.1100000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3g-'
      cbnlocal = '    1    7'
      canlocal = '    1    7'
      ctnlocal = '    1    7'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.9090000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4110000E+02
      nbclocal(5) =  -0.1620000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3h-'
      cbnlocal = '    1    8'
      canlocal = '    1    8'
      ctnlocal = '    1    8'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.9090000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4110000E+02
      nbclocal(5) =  -0.1620000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3i-'
      cbnlocal = '    1    9'
      canlocal = '    1    9'
      ctnlocal = '    1    9'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.1221000E+03
      nbclocal(3) =  -0.8320000E+02
      nbclocal(4) =  -0.4420000E+02
      nbclocal(5) =  -0.5300000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3j-'
      cbnlocal = '    1   10'
      canlocal = '    1   10'
      ctnlocal = '    1   10'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.1280000E+03
      nbclocal(3) =  -0.8910000E+02
      nbclocal(4) =  -0.5030000E+02
      nbclocal(5) =  -0.1140000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3k-'
      cbnlocal = '    1   11'
      canlocal = '    1   11'
      ctnlocal = '    1   11'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.9090000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4110000E+02
      nbclocal(5) =  -0.1620000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3l-'
      cbnlocal = '    1   12'
      canlocal = '    1   12'
      ctnlocal = '    1   12'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.4700000E+02
      nbclocal(3) =  -0.4380000E+02
      nbclocal(4) =  -0.4050000E+02
      nbclocal(5) =  -0.3730000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3m-'
      cbnlocal = '    1   13'
      canlocal = '    1   13'
      ctnlocal = '    1   13'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.9100000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4100000E+02
      nbclocal(5) =  -0.1600000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3n-'
      cbnlocal = '    1   14'
      canlocal = '    1   14'
      ctnlocal = '    1   14'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.3500000E+02
      nbclocal(3) =  -0.3470000E+02
      nbclocal(4) =  -0.3430000E+02
      nbclocal(5) =  -0.3400000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3o-'
      cbnlocal = '    1   15'
      canlocal = '    1   15'
      ctnlocal = '    1   15'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.3800000E+02
      nbclocal(3) =  -0.3510000E+02
      nbclocal(4) =  -0.3210000E+02
      nbclocal(5) =  -0.2920000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3p-'
      cbnlocal = '    1   16'
      canlocal = '    1   16'
      ctnlocal = '    1   16'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.4700000E+02
      nbclocal(3) =  -0.4700000E+02
      nbclocal(4) =  -0.4700000E+02
      nbclocal(5) =  -0.4700000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3q-'
      cbnlocal = '    1   17'
      canlocal = '    1   17'
      ctnlocal = '    1   17'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5800000E+02
      nbclocal(4) =  -0.4600000E+02
      nbclocal(5) =  -0.3400000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3r-'
      cbnlocal = '    1   18'
      canlocal = '    1   18'
      ctnlocal = '    1   18'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.6950000E+02
      nbclocal(3) =  -0.5870000E+02
      nbclocal(4) =  -0.4780000E+02
      nbclocal(5) =  -0.3700000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3s-'
      cbnlocal = '    1   19'
      canlocal = '    1   19'
      ctnlocal = '    1   19'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.9090000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4110000E+02
      nbclocal(5) =  -0.1620000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3t-'
      cbnlocal = '    1   20'
      canlocal = '    1   20'
      ctnlocal = '    1   20'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.9090000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4110000E+02
      nbclocal(5) =  -0.1620000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3u-'
      cbnlocal = '    1   21'
      canlocal = '    1   21'
      ctnlocal = '    1   21'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.9090000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4110000E+02
      nbclocal(5) =  -0.1620000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.03490
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3v-'
      cbnlocal = '    1   22'
      canlocal = '    1   22'
      ctnlocal = '    1   22'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.1140000E+03
      nbclocal(3) =  -0.7820000E+02
      nbclocal(4) =  -0.4240000E+02
      nbclocal(5) =  -0.6550000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2    1'
      canlocal = '    2    1'
      ctnlocal = '    2    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2490000E+02
      nbclocal(4) =  -0.2420000E+02
      nbclocal(5) =  -0.2350000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2    2'
      canlocal = '    2    2'
      ctnlocal = '    2    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.3290000E+02
      nbclocal(3) =  -0.2850000E+02
      nbclocal(4) =  -0.2400000E+02
      nbclocal(5) =  -0.1960000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2    3'
      canlocal = '    2    3'
      ctnlocal = '    2    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2490000E+02
      nbclocal(4) =  -0.2420000E+02
      nbclocal(5) =  -0.2350000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2    4'
      canlocal = '    2    4'
      ctnlocal = '    2    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2490000E+02
      nbclocal(4) =  -0.2420000E+02
      nbclocal(5) =  -0.2350000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2    5'
      canlocal = '    2    5'
      ctnlocal = '    2    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.3310000E+02
      nbclocal(3) =  -0.3030000E+02
      nbclocal(4) =  -0.2740000E+02
      nbclocal(5) =  -0.2460000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2    6'
      canlocal = '    2    6'
      ctnlocal = '    2    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3620000E+01
      nbclocal(2) =  -0.2200000E+02
      nbclocal(3) =  -0.2100000E+02
      nbclocal(4) =  -0.2000000E+02
      nbclocal(5) =  -0.1900000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2    7'
      canlocal = '    2    7'
      ctnlocal = '    2    7'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2490000E+02
      nbclocal(4) =  -0.2420000E+02
      nbclocal(5) =  -0.2350000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2    8'
      canlocal = '    2    8'
      ctnlocal = '    2    8'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2500000E+02
      nbclocal(3) =  -0.2300000E+02
      nbclocal(4) =  -0.2100000E+02
      nbclocal(5) =  -0.1900000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2    9'
      canlocal = '    2    9'
      ctnlocal = '    2    9'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.3000000E+02
      nbclocal(3) =  -0.2700000E+02
      nbclocal(4) =  -0.2400000E+02
      nbclocal(5) =  -0.2100000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   10'
      canlocal = '    2   10'
      ctnlocal = '    2   10'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.4680000E+02
      nbclocal(3) =  -0.4210000E+02
      nbclocal(4) =  -0.3750000E+02
      nbclocal(5) =  -0.3280000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   11'
      canlocal = '    2   11'
      ctnlocal = '    2   11'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.3600000E+02
      nbclocal(3) =  -0.2870000E+02
      nbclocal(4) =  -0.2130000E+02
      nbclocal(5) =  -0.1400000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   12'
      canlocal = '    2   12'
      ctnlocal = '    2   12'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.3310000E+02
      nbclocal(3) =  -0.3030000E+02
      nbclocal(4) =  -0.2740000E+02
      nbclocal(5) =  -0.2460000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   13'
      canlocal = '    2   13'
      ctnlocal = '    2   13'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2500000E+02
      nbclocal(3) =  -0.2470000E+02
      nbclocal(4) =  -0.2430000E+02
      nbclocal(5) =  -0.2400000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   14'
      canlocal = '    2   14'
      ctnlocal = '    2   14'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.3310000E+02
      nbclocal(3) =  -0.3030000E+02
      nbclocal(4) =  -0.2740000E+02
      nbclocal(5) =  -0.2460000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   15'
      canlocal = '    2   15'
      ctnlocal = '    2   15'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2100000E+02
      nbclocal(3) =  -0.2100000E+02
      nbclocal(4) =  -0.2100000E+02
      nbclocal(5) =  -0.2100000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   16'
      canlocal = '    2   16'
      ctnlocal = '    2   16'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.3200000E+02
      nbclocal(3) =  -0.3200000E+02
      nbclocal(4) =  -0.3200000E+02
      nbclocal(5) =  -0.3200000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   17'
      canlocal = '    2   17'
      ctnlocal = '    2   17'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3620000E+01
      nbclocal(2) =  -0.5500000E+02
      nbclocal(3) =  -0.3270000E+02
      nbclocal(4) =  -0.3230000E+02
      nbclocal(5) =  -0.3500000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   18'
      canlocal = '    2   18'
      ctnlocal = '    2   18'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2490000E+02
      nbclocal(4) =  -0.2420000E+02
      nbclocal(5) =  -0.2350000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   19'
      canlocal = '    2   19'
      ctnlocal = '    2   19'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.3000000E+02
      nbclocal(3) =  -0.2700000E+02
      nbclocal(4) =  -0.2400000E+02
      nbclocal(5) =  -0.2100000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   20'
      canlocal = '    2   20'
      ctnlocal = '    2   20'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2490000E+02
      nbclocal(4) =  -0.2420000E+02
      nbclocal(5) =  -0.2350000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   21'
      canlocal = '    2   21'
      ctnlocal = '    2   21'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2490000E+02
      nbclocal(4) =  -0.2420000E+02
      nbclocal(5) =  -0.2350000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02690
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH2 '
      cbnlocal = '    2   22'
      canlocal = '    2   22'
      ctnlocal = '    2   22'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3570000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2490000E+02
      nbclocal(4) =  -0.2420000E+02
      nbclocal(5) =  -0.2350000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3    1'
      canlocal = '    3    1'
      ctnlocal = '    3    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.6800000E+01
      nbclocal(3) =  -0.6400000E+01
      nbclocal(4) =  -0.6100000E+01
      nbclocal(5) =  -0.5700000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3    2'
      canlocal = '    3    2'
      ctnlocal = '    3    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.1330000E+02
      nbclocal(3) =  -0.1190000E+02
      nbclocal(4) =  -0.1050000E+02
      nbclocal(5) =  -0.9100000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3    3'
      canlocal = '    3    3'
      ctnlocal = '    3    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3    4'
      canlocal = '    3    4'
      ctnlocal = '    3    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3    5'
      canlocal = '    3    5'
      ctnlocal = '    3    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.4030000E+02
      nbclocal(3) =  -0.3620000E+02
      nbclocal(4) =  -0.3210000E+02
      nbclocal(5) =  -0.2800000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3    6'
      canlocal = '    3    6'
      ctnlocal = '    3    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.1800000E+02
      nbclocal(3) =  -0.1420000E+02
      nbclocal(4) =  -0.1030000E+02
      nbclocal(5) =  -0.6500000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3    7'
      canlocal = '    3    7'
      ctnlocal = '    3    7'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3    8'
      canlocal = '    3    8'
      ctnlocal = '    3    8'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.2200000E+02
      nbclocal(3) =  -0.1870000E+02
      nbclocal(4) =  -0.1530000E+02
      nbclocal(5) =  -0.1200000E+02
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3    9'
      canlocal = '    3    9'
      ctnlocal = '    3    9'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   10'
      canlocal = '    3   10'
      ctnlocal = '    3   10'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   11'
      canlocal = '    3   11'
      ctnlocal = '    3   11'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   12'
      canlocal = '    3   12'
      ctnlocal = '    3   12'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   13'
      canlocal = '    3   13'
      ctnlocal = '    3   13'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   14'
      canlocal = '    3   14'
      ctnlocal = '    3   14'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   15'
      canlocal = '    3   15'
      ctnlocal = '    3   15'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   16'
      canlocal = '    3   16'
      ctnlocal = '    3   16'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   17'
      canlocal = '    3   17'
      ctnlocal = '    3   17'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   18'
      canlocal = '    3   18'
      ctnlocal = '    3   18'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   19'
      canlocal = '    3   19'
      ctnlocal = '    3   19'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>CH- '
      cbnlocal = '    3   20'
      canlocal = '    3   20'
      ctnlocal = '    3   20'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3900000E+01
      nbclocal(2) =  -0.3150000E+02
      nbclocal(3) =  -0.2250000E+02
      nbclocal(4) =  -0.1340000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>C<  '
      cbnlocal = '    4    1'
      canlocal = '    4    1'
      ctnlocal = '    4    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.8500000E+01
      nbclocal(3) =  -0.8000000E+01
      nbclocal(4) =  -0.7600000E+01
      nbclocal(5) =  -0.7100000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>C<  '
      cbnlocal = '    4    2'
      canlocal = '    4    2'
      ctnlocal = '    4    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.8500000E+01
      nbclocal(3) =  -0.8000000E+01
      nbclocal(4) =  -0.7600000E+01
      nbclocal(5) =  -0.7100000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>C<  '
      cbnlocal = '    4    3'
      canlocal = '    4    3'
      ctnlocal = '    4    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.8500000E+01
      nbclocal(3) =  -0.8000000E+01
      nbclocal(4) =  -0.7600000E+01
      nbclocal(5) =  -0.7100000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>C<  '
      cbnlocal = '    4    4'
      canlocal = '    4    4'
      ctnlocal = '    4    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.8500000E+01
      nbclocal(3) =  -0.8000000E+01
      nbclocal(4) =  -0.7600000E+01
      nbclocal(5) =  -0.7100000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>C<  '
      cbnlocal = '    4    5'
      canlocal = '    4    5'
      ctnlocal = '    4    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.8500000E+01
      nbclocal(3) =  -0.8000000E+01
      nbclocal(4) =  -0.7600000E+01
      nbclocal(5) =  -0.7100000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>C<  '
      cbnlocal = '    4    6'
      canlocal = '    4    6'
      ctnlocal = '    4    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.8000000E+01
      nbclocal(3) =  -0.7700000E+01
      nbclocal(4) =  -0.7300000E+01
      nbclocal(5) =  -0.7000000E+01
      nbclocal(6) =   0.7700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'ACH  '
      cbnlocal = '    5    1'
      canlocal = '    5    1'
      ctnlocal = '    5    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.3100000E+02
      nbclocal(3) =  -0.2790000E+02
      nbclocal(4) =  -0.2480000E+02
      nbclocal(5) =  -0.2170000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'ACH  '
      cbnlocal = '    5    2'
      canlocal = '    5    2'
      ctnlocal = '    5    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.5500000E+02
      nbclocal(3) =  -0.4270000E+02
      nbclocal(4) =  -0.3030000E+02
      nbclocal(5) =  -0.1800000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'ACH  '
      cbnlocal = '    5    3'
      canlocal = '    5    3'
      ctnlocal = '    5    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.3100000E+02
      nbclocal(3) =  -0.2730000E+02
      nbclocal(4) =  -0.2370000E+02
      nbclocal(5) =  -0.1100000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'ACH  '
      cbnlocal = '    5    4'
      canlocal = '    5    4'
      ctnlocal = '    5    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.3100000E+02
      nbclocal(3) =  -0.3210000E+02
      nbclocal(4) =  -0.2460000E+02
      nbclocal(5) =  -0.2300000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'ACH  '
      cbnlocal = '    5    5'
      canlocal = '    5    5'
      ctnlocal = '    5    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3425000E+01
      nbclocal(2) =  -0.3100000E+02
      nbclocal(3) =  -0.3210000E+02
      nbclocal(4) =  -0.2460000E+02
      nbclocal(5) =  -0.2170000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'AC-  '
      cbnlocal = '    6    1'
      canlocal = '    6    1'
      ctnlocal = '    6    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3300000E+01
      nbclocal(2) =  -0.1700000E+02
      nbclocal(3) =  -0.1670000E+02
      nbclocal(4) =  -0.1630000E+02
      nbclocal(5) =  -0.1600000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'AC-  '
      cbnlocal = '    6    2'
      canlocal = '    6    2'
      ctnlocal = '    6    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3300000E+01
      nbclocal(2) =  -0.1500000E+02
      nbclocal(3) =  -0.1500000E+02
      nbclocal(4) =  -0.1500000E+02
      nbclocal(5) =  -0.1500000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'AC-  '
      cbnlocal = '    6    3'
      canlocal = '    6    3'
      ctnlocal = '    6    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3300000E+01
      nbclocal(2) =  -0.3500000E+02
      nbclocal(3) =  -0.2770000E+02
      nbclocal(4) =  -0.2030000E+02
      nbclocal(5) =  -0.1300000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'AC-  '
      cbnlocal = '    6    4'
      canlocal = '    6    4'
      ctnlocal = '    6    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3300000E+01
      nbclocal(2) =  -0.3500000E+02
      nbclocal(3) =  -0.2770000E+02
      nbclocal(4) =  -0.2030000E+02
      nbclocal(5) =  -0.1300000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'AC-  '
      cbnlocal = '    6    5'
      canlocal = '    6    5'
      ctnlocal = '    6    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3300000E+01
      nbclocal(2) =  -0.2500000E+02
      nbclocal(3) =  -0.2130000E+02
      nbclocal(4) =  -0.1770000E+02
      nbclocal(5) =  -0.1400000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02700
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2=a'
      cbnlocal = '    7    1'
      canlocal = '    7    1'
      ctnlocal = '    7    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.7570000E+02
      nbclocal(3) =  -0.5479000E+02
      nbclocal(4) =  -0.3389000E+02
      nbclocal(5) =  -0.1300000E+02
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02700
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2=b'
      cbnlocal = '    7    2'
      canlocal = '    7    2'
      ctnlocal = '    7    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.9090000E+02
      nbclocal(3) =  -0.6300000E+02
      nbclocal(4) =  -0.3500000E+02
      nbclocal(5) =  -0.7100000E+01
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02700
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2=c'
      cbnlocal = '    7    3'
      canlocal = '    7    3'
      ctnlocal = '    7    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.6500000E+02
      nbclocal(3) =  -0.5100000E+02
      nbclocal(4) =  -0.3700000E+02
      nbclocal(5) =  -0.2300000E+02
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02700
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2=b'
      cbnlocal = '    7    4'
      canlocal = '    7    4'
      ctnlocal = '    7    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.6500000E+02
      nbclocal(3) =  -0.5100000E+02
      nbclocal(4) =  -0.3700000E+02
      nbclocal(5) =  -0.2300000E+02
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.02700
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2=b'
      cbnlocal = '    7    5'
      canlocal = '    7    5'
      ctnlocal = '    7    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.6500000E+02
      nbclocal(3) =  -0.5100000E+02
      nbclocal(4) =  -0.3700000E+02
      nbclocal(5) =  -0.2300000E+02
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=CH- '
      cbnlocal = '    8    1'
      canlocal = '    8    1'
      ctnlocal = '    8    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2430000E+02
      nbclocal(4) =  -0.2290000E+02
      nbclocal(5) =  -0.2160000E+02
      nbclocal(6) =   0.7050000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=CH- '
      cbnlocal = '    8    2'
      canlocal = '    8    2'
      ctnlocal = '    8    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.4000000E+02
      nbclocal(3) =  -0.3320000E+02
      nbclocal(4) =  -0.2630000E+02
      nbclocal(5) =  -0.1950000E+02
      nbclocal(6) =   0.7050000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=CH- '
      cbnlocal = '    8    3'
      canlocal = '    8    3'
      ctnlocal = '    8    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2390000E+02
      nbclocal(4) =  -0.2220000E+02
      nbclocal(5) =  -0.2050000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=CH- '
      cbnlocal = '    8    4'
      canlocal = '    8    4'
      ctnlocal = '    8    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2390000E+02
      nbclocal(4) =  -0.2220000E+02
      nbclocal(5) =  -0.2050000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=CH- '
      cbnlocal = '    8    5'
      canlocal = '    8    5'
      ctnlocal = '    8    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2390000E+02
      nbclocal(4) =  -0.2220000E+02
      nbclocal(5) =  -0.2250000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=CH- '
      cbnlocal = '    8    6'
      canlocal = '    8    6'
      ctnlocal = '    8    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.2560000E+02
      nbclocal(3) =  -0.2390000E+02
      nbclocal(4) =  -0.2220000E+02
      nbclocal(5) =  -0.2250000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C<  '
      cbnlocal = '    9    1'
      canlocal = '    9    1'
      ctnlocal = '    9    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.7320000E+02
      nbclocal(3) =  -0.5660000E+02
      nbclocal(4) =  -0.3990000E+02
      nbclocal(5) =  -0.2330000E+02
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C<  '
      cbnlocal = '    9    2'
      canlocal = '    9    2'
      ctnlocal = '    9    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.6100000E+01
      nbclocal(3) =  -0.5200000E+01
      nbclocal(4) =  -0.4200000E+01
      nbclocal(5) =  -0.3300000E+01
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C<  '
      cbnlocal = '    9    3'
      canlocal = '    9    3'
      ctnlocal = '    9    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.1400000E+02
      nbclocal(3) =  -0.1200000E+02
      nbclocal(4) =  -0.1000000E+02
      nbclocal(5) =  -0.8000000E+01
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C<  '
      cbnlocal = '    9    4'
      canlocal = '    9    4'
      ctnlocal = '    9    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.1080000E+02
      nbclocal(3) =  -0.7400000E+01
      nbclocal(4) =  -0.4000000E+01
      nbclocal(5) =  -0.6000000E+00
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C<  '
      cbnlocal = '    9    5'
      canlocal = '    9    5'
      ctnlocal = '    9    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.1080000E+02
      nbclocal(3) =  -0.7400000E+01
      nbclocal(4) =  -0.4000000E+01
      nbclocal(5) =  -0.6000000E+00
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C=  '
      cbnlocal = '   10    1'
      canlocal = '   10    1'
      ctnlocal = '   10    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.3300000E+02
      nbclocal(3) =  -0.2440000E+02
      nbclocal(4) =  -0.1580000E+02
      nbclocal(5) =  -0.7200000E+01
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C=  '
      cbnlocal = '   10    2'
      canlocal = '   10    2'
      ctnlocal = '   10    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.3300000E+02
      nbclocal(3) =  -0.2440000E+02
      nbclocal(4) =  -0.1580000E+02
      nbclocal(5) =  -0.7200000E+01
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C=  '
      cbnlocal = '   10    3'
      canlocal = '   10    3'
      ctnlocal = '   10    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.3300000E+02
      nbclocal(3) =  -0.2440000E+02
      nbclocal(4) =  -0.1580000E+02
      nbclocal(5) =  -0.7200000E+01
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C=  '
      cbnlocal = '   10    4'
      canlocal = '   10    4'
      ctnlocal = '   10    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.3300000E+02
      nbclocal(3) =  -0.2440000E+02
      nbclocal(4) =  -0.1580000E+02
      nbclocal(5) =  -0.7200000E+01
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=C=  '
      cbnlocal = '   10    5'
      canlocal = '   10    5'
      ctnlocal = '   10    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3450000E+01
      nbclocal(2) =  -0.3300000E+02
      nbclocal(3) =  -0.2440000E+02
      nbclocal(4) =  -0.1580000E+02
      nbclocal(5) =  -0.7200000E+01
      nbclocal(6) =   0.6700000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<CH  '
      cbnlocal = '   11    1'
      canlocal = '   11    1'
      ctnlocal = '   11    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3300000E+01
      nbclocal(2) =  -0.9991000E+02
      nbclocal(3) =  -0.6780000E+02
      nbclocal(4) =  -0.3560000E+02
      nbclocal(5) =  -0.3500000E+01
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<CH  '
      cbnlocal = '   11    2'
      canlocal = '   11    2'
      ctnlocal = '   11    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.1156000E+03
      nbclocal(3) =  -0.8220000E+02
      nbclocal(4) =  -0.4880000E+02
      nbclocal(5) =  -0.1540000E+02
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<CH  '
      cbnlocal = '   11    3'
      canlocal = '   11    3'
      ctnlocal = '   11    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.1156000E+03
      nbclocal(3) =  -0.8220000E+02
      nbclocal(4) =  -0.4880000E+02
      nbclocal(5) =  -0.1540000E+02
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<CH  '
      cbnlocal = '   11    4'
      canlocal = '   11    4'
      ctnlocal = '   11    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.1156000E+03
      nbclocal(3) =  -0.8220000E+02
      nbclocal(4) =  -0.4880000E+02
      nbclocal(5) =  -0.1540000E+02
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   13.01910
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<CH  '
      cbnlocal = '   11    5'
      canlocal = '   11    5'
      ctnlocal = '   11    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.1156000E+03
      nbclocal(3) =  -0.8220000E+02
      nbclocal(4) =  -0.4880000E+02
      nbclocal(5) =  -0.1540000E+02
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<C-  '
      cbnlocal = '   12    1'
      canlocal = '   12    1'
      ctnlocal = '   12    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.5145000E+02
      nbclocal(3) =  -0.3940000E+02
      nbclocal(4) =  -0.2740000E+02
      nbclocal(5) =  -0.1530000E+02
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<C-  '
      cbnlocal = '   12    2'
      canlocal = '   12    2'
      ctnlocal = '   12    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.1050000E+03
      nbclocal(3) =  -0.8700000E+02
      nbclocal(4) =  -0.6900000E+02
      nbclocal(5) =  -0.5100000E+02
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<C-  '
      cbnlocal = '   12    3'
      canlocal = '   12    3'
      ctnlocal = '   12    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3500000E+01
      nbclocal(2) =  -0.8150000E+02
      nbclocal(3) =  -0.5630000E+02
      nbclocal(4) =  -0.3100000E+02
      nbclocal(5) =  -0.5800000E+01
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<C-  '
      cbnlocal = '   12    4'
      canlocal = '   12    4'
      ctnlocal = '   12    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3200000E+01
      nbclocal(2) =  -0.7200000E+02
      nbclocal(3) =  -0.5470000E+02
      nbclocal(4) =  -0.3730000E+02
      nbclocal(5) =  -0.2000000E+02
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   12.01120
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<C-  '
      cbnlocal = '   12    5'
      canlocal = '   12    5'
      ctnlocal = '   12    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3200000E+01
      nbclocal(2) =  -0.7200000E+02
      nbclocal(3) =  -0.5470000E+02
      nbclocal(4) =  -0.3730000E+02
      nbclocal(5) =  -0.2000000E+02
      nbclocal(6) =   0.6000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   18.01520
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2O  '
      cbnlocal = '   13    1'
      canlocal = '   13    1'
      ctnlocal = '   13    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3000000E+01
      nbclocal(2) =  -0.2610000E+03
      nbclocal(3) =  -0.1755000E+03
      nbclocal(4) =  -0.8990000E+02
      nbclocal(5) =  -0.4400000E+01
      nbclocal(6) =   0.6300000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '(Me)O'
      cbnlocal = '   14    1'
      canlocal = '   14    1'
      ctnlocal = '   14    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1725000E+03
      nbclocal(3) =  -0.1238000E+03
      nbclocal(4) =  -0.7510000E+02
      nbclocal(5) =  -0.2640000E+02
      nbclocal(6) =   0.6300000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14    2'
      canlocal = '   14    2'
      ctnlocal = '   14    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1382000E+03
      nbclocal(3) =  -0.9690000E+02
      nbclocal(4) =  -0.5550000E+02
      nbclocal(5) =  -0.1420000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14    3'
      canlocal = '   14    3'
      ctnlocal = '   14    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1715000E+03
      nbclocal(3) =  -0.1180000E+03
      nbclocal(4) =  -0.6440000E+02
      nbclocal(5) =  -0.1090000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14    4'
      canlocal = '   14    4'
      ctnlocal = '   14    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1248000E+03
      nbclocal(3) =  -0.1151000E+03
      nbclocal(4) =  -0.1053000E+03
      nbclocal(5) =  -0.9560000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14    5'
      canlocal = '   14    5'
      ctnlocal = '   14    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1110000E+03
      nbclocal(3) =  -0.7570000E+02
      nbclocal(4) =  -0.4040000E+02
      nbclocal(5) =  -0.5100000E+01
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14    6'
      canlocal = '   14    6'
      ctnlocal = '   14    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.9800000E+02
      nbclocal(3) =  -0.9730000E+02
      nbclocal(4) =  -0.9670000E+02
      nbclocal(5) =  -0.9600000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14    7'
      canlocal = '   14    7'
      ctnlocal = '   14    7'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1110000E+03
      nbclocal(3) =  -0.7570000E+02
      nbclocal(4) =  -0.4040000E+02
      nbclocal(5) =  -0.5100000E+01
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14    8'
      canlocal = '   14    8'
      ctnlocal = '   14    8'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1110000E+03
      nbclocal(3) =  -0.7570000E+02
      nbclocal(4) =  -0.4040000E+02
      nbclocal(5) =  -0.5100000E+01
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14    9'
      canlocal = '   14    9'
      ctnlocal = '   14    9'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1110000E+03
      nbclocal(3) =  -0.7570000E+02
      nbclocal(4) =  -0.4040000E+02
      nbclocal(5) =  -0.5100000E+01
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14   10'
      canlocal = '   14   10'
      ctnlocal = '   14   10'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1850000E+03
      nbclocal(3) =  -0.1262000E+03
      nbclocal(4) =  -0.6740000E+02
      nbclocal(5) =  -0.8600000E+01
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14   11'
      canlocal = '   14   11'
      ctnlocal = '   14   11'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2850000E+01
      nbclocal(2) =  -0.1850000E+03
      nbclocal(3) =  -0.1262000E+03
      nbclocal(4) =  -0.6740000E+02
      nbclocal(5) =  -0.8600000E+01
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   17.00730
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH   '
      cbnlocal = '   14   12'
      canlocal = '   14   12'
      ctnlocal = '   14   12'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2600000E+01
      nbclocal(2) =  -0.6100000E+02
      nbclocal(3) =  -0.1262000E+03
      nbclocal(4) =  -0.6740000E+02
      nbclocal(5) =  -0.5500000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15    1'
      canlocal = '   15    1'
      ctnlocal = '   15    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.5200000E+02
      nbclocal(3) =  -0.4770000E+02
      nbclocal(4) =  -0.4330000E+02
      nbclocal(5) =  -0.3900000E+02
      nbclocal(6) =   0.6800000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15    2'
      canlocal = '   15    2'
      ctnlocal = '   15    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.9120000E+02
      nbclocal(3) =  -0.6190000E+02
      nbclocal(4) =  -0.3250000E+02
      nbclocal(5) =  -0.3200000E+01
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15    3'
      canlocal = '   15    3'
      ctnlocal = '   15    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15    4'
      canlocal = '   15    4'
      ctnlocal = '   15    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.2874000E+03
      nbclocal(3) =  -0.2005000E+03
      nbclocal(4) =  -0.1136000E+03
      nbclocal(5) =  -0.2670000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15    5'
      canlocal = '   15    5'
      ctnlocal = '   15    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15    6'
      canlocal = '   15    6'
      ctnlocal = '   15    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15    7'
      canlocal = '   15    7'
      ctnlocal = '   15    7'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.2874000E+03
      nbclocal(3) =  -0.2005000E+03
      nbclocal(4) =  -0.1136000E+03
      nbclocal(5) =  -0.2670000E+02
      nbclocal(6) =   0.6600000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15    8'
      canlocal = '   15    8'
      ctnlocal = '   15    8'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.7800000E+02
      nbclocal(3) =  -0.6270000E+02
      nbclocal(4) =  -0.4730000E+02
      nbclocal(5) =  -0.3200000E+02
      nbclocal(6) =   0.6500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15    9'
      canlocal = '   15    9'
      ctnlocal = '   15    9'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O<   '
      cbnlocal = '   15   10'
      canlocal = '   15   10'
      ctnlocal = '   15   10'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16    1'
      canlocal = '   16    1'
      ctnlocal = '   16    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3120000E+01
      nbclocal(2) =  -0.6900000E+02
      nbclocal(3) =  -0.6630000E+02
      nbclocal(4) =  -0.6370000E+02
      nbclocal(5) =  -0.6100000E+02
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16    2'
      canlocal = '   16    2'
      ctnlocal = '   16    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.1490000E+03
      nbclocal(3) =  -0.1218000E+03
      nbclocal(4) =  -0.9460000E+02
      nbclocal(5) =  -0.6740000E+02
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16    3'
      canlocal = '   16    3'
      ctnlocal = '   16    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2600000E+01
      nbclocal(2) =  -0.6500000E+02
      nbclocal(3) =  -0.1218000E+03
      nbclocal(4) =  -0.9460000E+02
      nbclocal(5) =  -0.5840000E+02
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16    4'
      canlocal = '   16    4'
      ctnlocal = '   16    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2900000E+01
      nbclocal(2) =  -0.4200000E+02
      nbclocal(3) =  -0.3830000E+02
      nbclocal(4) =  -0.3470000E+02
      nbclocal(5) =  -0.3100000E+02
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16    5'
      canlocal = '   16    5'
      ctnlocal = '   16    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2900000E+01
      nbclocal(2) =  -0.1410000E+03
      nbclocal(3) =  -0.1310000E+03
      nbclocal(4) =  -0.1210000E+03
      nbclocal(5) =  -0.1110000E+03
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16    6'
      canlocal = '   16    6'
      ctnlocal = '   16    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3120000E+01
      nbclocal(2) =  -0.9700000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.5800000E+02
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16    7'
      canlocal = '   16    7'
      ctnlocal = '   16    7'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16    8'
      canlocal = '   16    8'
      ctnlocal = '   16    8'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16    9'
      canlocal = '   16    9'
      ctnlocal = '   16    9'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.99940
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=   '
      cbnlocal = '   16   10'
      canlocal = '   16   10'
      ctnlocal = '   16   10'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2700000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6200000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   32.06620
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S<   '
      cbnlocal = '   17    1'
      canlocal = '   17    1'
      ctnlocal = '   17    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3260000E+01
      nbclocal(2) =  -0.1110000E+03
      nbclocal(3) =  -0.1037000E+03
      nbclocal(4) =  -0.9630000E+02
      nbclocal(5) =  -0.8900000E+02
      nbclocal(6) =   0.1040000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   32.06620
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S<   '
      cbnlocal = '   17    2'
      canlocal = '   17    2'
      ctnlocal = '   17    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3260000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1040000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   32.06620
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S<   '
      cbnlocal = '   17    3'
      canlocal = '   17    3'
      ctnlocal = '   17    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3260000E+01
      nbclocal(2) =  -0.1110000E+03
      nbclocal(3) =  -0.9770000E+02
      nbclocal(4) =  -0.8430000E+02
      nbclocal(5) =  -0.7000000E+02
      nbclocal(6) =   0.1040000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   32.06620
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S<   '
      cbnlocal = '   17    4'
      canlocal = '   17    4'
      ctnlocal = '   17    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3850000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1040000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   38.07180
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S<   '
      cbnlocal = '   17    5'
      canlocal = '   17    5'
      ctnlocal = '   17    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3630000E+01
      nbclocal(2) =  -0.1900000E+03
      nbclocal(3) =  -0.1680000E+03
      nbclocal(4) =  -0.9500000E+02
      nbclocal(5) =  -0.4440000E+02
      nbclocal(6) =   0.1160000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   33.07480
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-SH  '
      cbnlocal = '   18    1'
      canlocal = '   18    1'
      ctnlocal = '   18    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3600000E+01
      nbclocal(2) =  -0.1088000E+03
      nbclocal(3) =  -0.1063000E+03
      nbclocal(4) =  -0.1038000E+03
      nbclocal(5) =  -0.1013000E+03
      nbclocal(6) =   0.1040000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   33.07480
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-SH  '
      cbnlocal = '   18    2'
      canlocal = '   18    2'
      ctnlocal = '   18    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.4000000E+01
      nbclocal(2) =  -0.1110000E+03
      nbclocal(3) =  -0.9300000E+02
      nbclocal(4) =  -0.7500000E+02
      nbclocal(5) =  -0.5700000E+02
      nbclocal(6) =   0.1040000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   33.07480
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-SH  '
      cbnlocal = '   18    3'
      canlocal = '   18    3'
      ctnlocal = '   18    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1040000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   33.07480
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-SH  '
      cbnlocal = '   18    4'
      canlocal = '   18    4'
      ctnlocal = '   18    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1040000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   33.07480
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-S   '
      cbnlocal = '   18    5'
      canlocal = '   18    5'
      ctnlocal = '   18    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1040000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   16.02270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH2 '
      cbnlocal = '   19    1'
      canlocal = '   19    1'
      ctnlocal = '   19    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.8500000E+02
      nbclocal(3) =  -0.8130000E+02
      nbclocal(4) =  -0.7770000E+02
      nbclocal(5) =  -0.7400000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   16.02270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH2 '
      cbnlocal = '   19    2'
      canlocal = '   19    2'
      ctnlocal = '   19    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   16.02270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH2 '
      cbnlocal = '   19    3'
      canlocal = '   19    3'
      ctnlocal = '   19    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.1410000E+03
      nbclocal(3) =  -0.1310000E+03
      nbclocal(4) =  -0.1210000E+03
      nbclocal(5) =  -0.1110000E+03
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   16.02270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH2 '
      cbnlocal = '   19    4'
      canlocal = '   19    4'
      ctnlocal = '   19    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   16.02270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH2 '
      cbnlocal = '   19    5'
      canlocal = '   19    5'
      ctnlocal = '   19    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.01470
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH- '
      cbnlocal = '   20    1'
      canlocal = '   20    1'
      ctnlocal = '   20    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.8500000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4700000E+02
      nbclocal(5) =  -0.2800000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.01470
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH- '
      cbnlocal = '   20    2'
      canlocal = '   20    2'
      ctnlocal = '   20    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.01470
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH- '
      cbnlocal = '   20    3'
      canlocal = '   20    3'
      ctnlocal = '   20    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.8500000E+02
      nbclocal(3) =  -0.6600000E+02
      nbclocal(4) =  -0.4700000E+02
      nbclocal(5) =  -0.2800000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.01470
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH- '
      cbnlocal = '   20    4'
      canlocal = '   20    4'
      ctnlocal = '   20    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   15.01470
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-NH- '
      cbnlocal = '   20    5'
      canlocal = '   20    5'
      ctnlocal = '   20    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-N<  '
      cbnlocal = '   21    1'
      canlocal = '   21    1'
      ctnlocal = '   21    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.4820000E+02
      nbclocal(4) =  -0.2630000E+02
      nbclocal(5) =  -0.4500000E+01
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-N<  '
      cbnlocal = '   21    2'
      canlocal = '   21    2'
      ctnlocal = '   21    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-N<  '
      cbnlocal = '   21    3'
      canlocal = '   21    3'
      ctnlocal = '   21    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.1080000E+02
      nbclocal(3) =  -0.1080000E+02
      nbclocal(4) =  -0.1080000E+02
      nbclocal(5) =  -0.1080000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-N<  '
      cbnlocal = '   21    4'
      canlocal = '   21    4'
      ctnlocal = '   21    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-N<  '
      cbnlocal = '   21    5'
      canlocal = '   21    5'
      ctnlocal = '   21    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.7000000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=AN- '
      cbnlocal = '   22    1'
      canlocal = '   22    1'
      ctnlocal = '   22    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=AN- '
      cbnlocal = '   22    2'
      canlocal = '   22    2'
      ctnlocal = '   22    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=AN- '
      cbnlocal = '   22    3'
      canlocal = '   22    3'
      ctnlocal = '   22    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=AN- '
      cbnlocal = '   22    4'
      canlocal = '   22    4'
      ctnlocal = '   22    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '=AN- '
      cbnlocal = '   22    5'
      canlocal = '   22    5'
      ctnlocal = '   22    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.6500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<N   '
      cbnlocal = '   23    1'
      canlocal = '   23    1'
      ctnlocal = '   23    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3080000E+01
      nbclocal(2) =  -0.3300000E+02
      nbclocal(3) =  -0.2450000E+02
      nbclocal(4) =  -0.1610000E+02
      nbclocal(5) =  -0.7600000E+01
      nbclocal(6) =   0.5500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '<N   '
      cbnlocal = '   23    2'
      canlocal = '   23    2'
      ctnlocal = '   23    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3080000E+01
      nbclocal(2) =  -0.1340000E+03
      nbclocal(3) =  -0.9120000E+02
      nbclocal(4) =  -0.4840000E+02
      nbclocal(5) =  -0.5600000E+01
      nbclocal(6) =   0.5500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '#N   '
      cbnlocal = '   23    3'
      canlocal = '   23    3'
      ctnlocal = '   23    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3080000E+01
      nbclocal(2) =  -0.3300000E+02
      nbclocal(3) =  -0.2450000E+02
      nbclocal(4) =  -0.1610000E+02
      nbclocal(5) =  -0.7600000E+01
      nbclocal(6) =   0.5500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '#N   '
      cbnlocal = '   23    4'
      canlocal = '   23    4'
      ctnlocal = '   23    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.5500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   14.00670
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '#N   '
      cbnlocal = '   23    5'
      canlocal = '   23    5'
      ctnlocal = '   23    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3100000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.5500000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   35.45270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Cl  '
      cbnlocal = '   24    1'
      canlocal = '   24    1'
      ctnlocal = '   24    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3330000E+01
      nbclocal(2) =  -0.1638000E+03
      nbclocal(3) =  -0.1147000E+03
      nbclocal(4) =  -0.6550000E+02
      nbclocal(5) =  -0.1640000E+02
      nbclocal(6) =   0.9940000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   35.45270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Cl  '
      cbnlocal = '   24    2'
      canlocal = '   24    2'
      ctnlocal = '   24    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3330000E+01
      nbclocal(2) =  -0.1740000E+03
      nbclocal(3) =  -0.9770000E+02
      nbclocal(4) =  -0.5230000E+02
      nbclocal(5) =  -0.3800000E+02
      nbclocal(6) =   0.9940000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   35.45270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Cl  '
      cbnlocal = '   24    3'
      canlocal = '   24    3'
      ctnlocal = '   24    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3200000E+01
      nbclocal(2) =  -0.1964000E+03
      nbclocal(3) =  -0.1375000E+03
      nbclocal(4) =  -0.7850000E+02
      nbclocal(5) =  -0.1960000E+02
      nbclocal(6) =   0.9940000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   35.45270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Cl  '
      cbnlocal = '   24    4'
      canlocal = '   24    4'
      ctnlocal = '   24    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3200000E+01
      nbclocal(2) =  -0.1638000E+03
      nbclocal(3) =  -0.1147000E+03
      nbclocal(4) =  -0.6550000E+02
      nbclocal(5) =  -0.1640000E+02
      nbclocal(6) =   0.9940000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   35.45270
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Cl  '
      cbnlocal = '   24    5'
      canlocal = '   24    5'
      ctnlocal = '   24    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3330000E+01
      nbclocal(2) =  -0.8430000E+02
      nbclocal(3) =  -0.5900000E+02
      nbclocal(4) =  -0.3370000E+02
      nbclocal(5) =  -0.8400000E+01
      nbclocal(6) =   0.9940000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   18.99840
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-F   '
      cbnlocal = '   25    1'
      canlocal = '   25    1'
      ctnlocal = '   25    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2640000E+01
      nbclocal(2) =  -0.2980000E+02
      nbclocal(3) =  -0.2690000E+02
      nbclocal(4) =  -0.2410000E+02
      nbclocal(5) =  -0.2120000E+02
      nbclocal(6) =   0.6400000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   18.99840
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-F   '
      cbnlocal = '   25    2'
      canlocal = '   25    2'
      ctnlocal = '   25    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2640000E+01
      nbclocal(2) =  -0.3180000E+02
      nbclocal(3) =  -0.2700000E+02
      nbclocal(4) =  -0.2220000E+02
      nbclocal(5) =  -0.1740000E+02
      nbclocal(6) =   0.6400000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   18.99840
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-F   '
      cbnlocal = '   25    3'
      canlocal = '   25    3'
      ctnlocal = '   25    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2640000E+01
      nbclocal(2) =  -0.2400000E+02
      nbclocal(3) =  -0.2700000E+02
      nbclocal(4) =  -0.2220000E+02
      nbclocal(5) =  -0.1740000E+02
      nbclocal(6) =   0.6400000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   18.99840
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-F   '
      cbnlocal = '   25    4'
      canlocal = '   25    4'
      ctnlocal = '   25    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2640000E+01
      nbclocal(2) =  -0.2400000E+02
      nbclocal(3) =  -0.2700000E+02
      nbclocal(4) =  -0.2220000E+02
      nbclocal(5) =  -0.1420000E+02
      nbclocal(6) =   0.6400000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   18.99840
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-F   '
      cbnlocal = '   25    5'
      canlocal = '   25    5'
      ctnlocal = '   25    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2640000E+01
      nbclocal(2) =  -0.3180000E+02
      nbclocal(3) =  -0.2700000E+02
      nbclocal(4) =  -0.2220000E+02
      nbclocal(5) =  -0.1740000E+02
      nbclocal(6) =   0.6400000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   18.99840
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-F   '
      cbnlocal = '   25    6'
      canlocal = '   25    6'
      ctnlocal = '   25    6'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.2640000E+01
      nbclocal(2) =  -0.3180000E+02
      nbclocal(3) =  -0.2700000E+02
      nbclocal(4) =  -0.2220000E+02
      nbclocal(5) =  -0.1740000E+02
      nbclocal(6) =   0.6400000E-01
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.97380
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>P   '
      cbnlocal = '   26    0'
      canlocal = '   26    0'
      ctnlocal = '   26    0'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.4000000E+01
      nbclocal(2) =  -0.2200000E+02
      nbclocal(3) =  -0.2170000E+02
      nbclocal(4) =  -0.2130000E+02
      nbclocal(5) =  -0.2100000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.97380
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>P<  '
      cbnlocal = '   26    1'
      canlocal = '   26    1'
      ctnlocal = '   26    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.4000000E+01
      nbclocal(2) =  -0.2200000E+02
      nbclocal(3) =  -0.2170000E+02
      nbclocal(4) =  -0.2130000E+02
      nbclocal(5) =  -0.2100000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.97380
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-P<  '
      cbnlocal = '   26    2'
      canlocal = '   26    2'
      ctnlocal = '   26    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.97380
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-P<  '
      cbnlocal = '   26    3'
      canlocal = '   26    3'
      ctnlocal = '   26    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.97380
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-P<  '
      cbnlocal = '   26    4'
      canlocal = '   26    4'
      ctnlocal = '   26    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.97380
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-P<  '
      cbnlocal = '   26    5'
      canlocal = '   26    5'
      ctnlocal = '   26    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   79.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Br  '
      cbnlocal = '   27    1'
      canlocal = '   27    1'
      ctnlocal = '   27    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3460000E+01
      nbclocal(2) =  -0.2272000E+03
      nbclocal(3) =  -0.1598000E+03
      nbclocal(4) =  -0.9250000E+02
      nbclocal(5) =  -0.2510000E+02
      nbclocal(6) =   0.1140000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   79.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Br  '
      cbnlocal = '   27    2'
      canlocal = '   27    2'
      ctnlocal = '   27    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   79.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Br  '
      cbnlocal = '   27    3'
      canlocal = '   27    3'
      ctnlocal = '   27    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   79.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Br  '
      cbnlocal = '   27    4'
      canlocal = '   27    4'
      ctnlocal = '   27    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   79.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-Br  '
      cbnlocal = '   27    5'
      canlocal = '   27    5'
      ctnlocal = '   27    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1100000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =  126.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-I   '
      cbnlocal = '   28    1'
      canlocal = '   28    1'
      ctnlocal = '   28    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1333000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =  126.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-I   '
      cbnlocal = '   28    2'
      canlocal = '   28    2'
      ctnlocal = '   28    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1800000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =  126.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-I   '
      cbnlocal = '   28    3'
      canlocal = '   28    3'
      ctnlocal = '   28    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1800000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =  126.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-I   '
      cbnlocal = '   28    4'
      canlocal = '   28    4'
      ctnlocal = '   28    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1800000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =  126.90400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-I   '
      cbnlocal = '   28    5'
      canlocal = '   28    5'
      ctnlocal = '   28    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1800000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   31.10990
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-SiH3'
      cbnlocal = '   29    1'
      canlocal = '   29    1'
      ctnlocal = '   29    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   31.10990
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-SiH3'
      cbnlocal = '   29    2'
      canlocal = '   29    2'
      ctnlocal = '   29    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   31.10990
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-SiH3'
      cbnlocal = '   29    3'
      canlocal = '   29    3'
      ctnlocal = '   29    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   31.10990
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-SiH3'
      cbnlocal = '   29    4'
      canlocal = '   29    4'
      ctnlocal = '   29    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   31.10990
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '-SiH3'
      cbnlocal = '   29    5'
      canlocal = '   29    5'
      ctnlocal = '   29    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.00610
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH2'
      cbnlocal = '   30    1'
      canlocal = '   30    1'
      ctnlocal = '   30    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.00610
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH2'
      cbnlocal = '   30    2'
      canlocal = '   30    2'
      ctnlocal = '   30    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.00610
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH2'
      cbnlocal = '   30    3'
      canlocal = '   30    3'
      ctnlocal = '   30    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.00610
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH2'
      cbnlocal = '   30    4'
      canlocal = '   30    4'
      ctnlocal = '   30    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   30.00610
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH2'
      cbnlocal = '   30    5'
      canlocal = '   30    5'
      ctnlocal = '   30    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   29.09400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH-'
      cbnlocal = '   31    1'
      canlocal = '   31    1'
      ctnlocal = '   31    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   29.09400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH-'
      cbnlocal = '   31    2'
      canlocal = '   31    2'
      ctnlocal = '   31    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   29.09400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH-'
      cbnlocal = '   31    3'
      canlocal = '   31    3'
      ctnlocal = '   31    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   29.09400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH-'
      cbnlocal = '   31    4'
      canlocal = '   31    4'
      ctnlocal = '   31    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   29.09400
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>SiH-'
      cbnlocal = '   31    5'
      canlocal = '   31    5'
      ctnlocal = '   31    5'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   28.08600
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>Si< '
      cbnlocal = '   32    1'
      canlocal = '   32    1'
      ctnlocal = '   32    1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   28.08600
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>Si< '
      cbnlocal = '   32    2'
      canlocal = '   32    2'
      ctnlocal = '   32    2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   28.08600
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>Si< '
      cbnlocal = '   32    3'
      canlocal = '   32    3'
      ctnlocal = '   32    3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   28.08600
      elelocal =  ''
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = '>Si< '
      cbnlocal = '   32    4'
      canlocal = '   32    4'
      ctnlocal = '   32    4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   0.3750000E+01
      nbclocal(2) =  -0.7000000E+02
      nbclocal(3) =  -0.5480000E+02
      nbclocal(4) =  -0.3970000E+02
      nbclocal(5) =  -0.2450000E+02
      nbclocal(6) =   0.1000000E+00
      call twh_set_nbcoeff(ii,nbclocal) 

      return
      end
