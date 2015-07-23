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
      subroutine twh_ffcharmm27x(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * contains some supplementary parameters that I created in order *
c     * to extend Charmm27 to some other systems                       *
c     * generally extrapolated from the official Charmm27 parameters,  *
c     * but just things I made up so I could do some more simulations  *
c     * so use at your own risk                                        *
c     *                                                                *
c     * originally written 02-22-2005 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer type,isame
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
      double precision onepi

      ffname = 'Charmm27'
      type = 0
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

c     --- no additional nonbond parameters at the moment

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Charmm27 bond potentials                                       *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.360d0
      vclocal(1) = (302.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN1 - CN5 [from CN1 - CN5G]
      isame = isame + 1
      vnlocal(isame,1) = 'CN1'
      vnlocal(isame,2) = 'CN5'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * charmm27 angles                                                *
c     ******************************************************************

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.50d0
      bclocal(1) = (140.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NN2 - CN5 - NN2 [from NN2B - CN5 - NN3G]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN2'
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
C     --- CN1 - CN5 - NN4 [from CN1 - CN5 - NN4]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'NN4'
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
C     --- CN5 - CN1 - NN2 [from CN5G - CN1 - NN2G]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'NN2'
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
C     --- CN5 - CN1 - ON1 [from CN5G - CN1 - ON1]
      isame = isame + 1
      bnlocal(isame,1) = 'CN5'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'ON1'
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
C     --- CN1 - CN5 - CN5 [from CN1 - CN5G - CN5]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'CN5'
      bnlocal(isame,3) = 'CN5'
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
C     --- CN1 - NN2 - CN1 [from CN1T - NN2U - CN1]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN1'
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
C     --- NN2 - CN1 - ON1 [from NN2G - CN1 - ON1]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'ON1'
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
C     --- NN2 - CN1 - NN2 [from NN2B - CN1T - NN2U]
      isame = isame + 1
      bnlocal(isame,1) = 'NN2'
      bnlocal(isame,2) = 'CN1'
      bnlocal(isame,3) = 'NN2'
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
C     --- CN1 - NN2 - CN5 [from CN1T - NN2U - CN1]
      isame = isame + 1
      bnlocal(isame,1) = 'CN1'
      bnlocal(isame,2) = 'NN2'
      bnlocal(isame,3) = 'CN5'
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
C     --- CT3 - CC - CT3 [taken from the CHARMM27 CT3-CT2-CT2]
c     --- published for acetone in Martin & Biddy Fluid Phase Equilib.
c     --- 236 53-57 (2005)
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CC'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return


c     ******************************************************************
c     * Charmm27 standard  dihedral                                    *
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)

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
c     --- X - CN1 - CN5 - X  [from X - CN1 - CN5G - X]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'CN1'
      tnlocal(isame,3) = 'CN5'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)


      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
