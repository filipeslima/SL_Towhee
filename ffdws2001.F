#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2007-2014 Marcus G. Martin                       *
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
      subroutine twh_ffdws2001(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the dws2001      *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 01-10-2007 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
c     --- character string scalars
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
c     --- character string arrays
      character*(FFNAMELEN) ffname
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
c     --- integer scalars
      integer ii,type,isame
c     --- double precision scalars
      double precision masslocal
c     --- double precision arrays
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)

      ffname = 'DWS2001'
      call twh_classical_potential(GLB_SET,'Drukker H-bond')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

c     References:
c     [jcp2001]
c     K. Drukker; G. Wu; G. C. Schatz;
c     "Model simulatinos of DNA denaturation dynamics",
c     J. Chem. Phys. 114, 579-590 (2001)
c     [pc2007]
c     personal communication bewteen G. Wu and M.G. Martin 01/2007

      type = 0
      
c     --- backbone sites [jcp2001]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.4d0
      nbclocal(2) = 150.0d0
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- no charges on this model
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 178.0d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal =  'BB'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'backbone'
      cbnlocal = 'backbone'
      canlocal = 'backbone'
      ctnlocal = 'backbone'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     --- A [jcp2001]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.0d0
      nbclocal(2) = 225.0d0
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- no charges on this model
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 135.0d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal =  'AA'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'A'
      cbnlocal = 'A'
      canlocal = 'A'
      ctnlocal = 'A'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     --- T [jcp2001]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.4d0
      nbclocal(2) = 225.0d0
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- no charges on this model
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 126.0d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal =  'TT'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'T'
      cbnlocal = 'T'
      canlocal = 'T'
      ctnlocal = 'T'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     --- C [jcp2001]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.4d0
      nbclocal(2) = 225.0d0
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- no charges on this model
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 111.0d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal =  'CC'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C'
      cbnlocal = 'C'
      canlocal = 'C'
      ctnlocal = 'C'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     --- G [jcp2001]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.0d0
      nbclocal(2) = 225.0d0
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- no charges on this model
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 151.0d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal =  'GG'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'G'
      cbnlocal = 'G'
      canlocal = 'G'
      ctnlocal = 'G'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     --- harmonic vibrations
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC_NONBOND)
      vclocal(0) = 3.05d0
      vclocal(1) = (100.0d0)*(0.5d0)*CON_KCALTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     backbone - A [jcp2001]
      isame = isame + 1
      vnlocal(isame,1) = 'backbone'
      vnlocal(isame,2) = 'A'
c     backbone - G [jcp2001]
      isame = isame + 1
      vnlocal(isame,1) = 'backbone'
      vnlocal(isame,2) = 'G'
c     --- set nvsame
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC_NONBOND)
      vclocal(0) = 2.70d0
      vclocal(1) = (100.0d0)*(0.5d0)*CON_KCALTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     backbone - T [jcp2001]
      isame = isame + 1
      vnlocal(isame,1) = 'backbone'
      vnlocal(isame,2) = 'T'
c     backbone - C [jcp2001]
      isame = isame + 1
      vnlocal(isame,1) = 'backbone'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC_NONBOND)
      vclocal(0) = 3.4d0
      vclocal(1) = (100.0d0)*(0.5d0)*CON_KCALTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     backbone - backbone [pc2007]
      isame = isame + 1
      vnlocal(isame,1) = 'backbone'
      vnlocal(isame,2) = 'backbone'
c     --- set nvsame
      call twh_set_vibnames(type,isame,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     --- harmonic cosine angles plus 1-3 nonbond
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS_OT)
      bclocal(0) = 90.0d0
      bclocal(1) = (80.0d0)*(0.5d0)*CON_KCALTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     backbone - backbone - A [jcp2001 and pc2007]
      isame = isame + 1
      bnlocal(isame,1) = 'backbone'
      bnlocal(isame,2) = 'backbone'
      bnlocal(isame,3) = 'A'
c     backbone - backbone - T [jcp2001 and pc2007]
      isame = isame + 1
      bnlocal(isame,1) = 'backbone'
      bnlocal(isame,2) = 'backbone'
      bnlocal(isame,3) = 'T'
c     backbone - backbone - C [jcp2001 and pc2007]
      isame = isame + 1
      bnlocal(isame,1) = 'backbone'
      bnlocal(isame,2) = 'backbone'
      bnlocal(isame,3) = 'C'
c     backbone - backbone - G [jcp2001 and pc2007]
      isame = isame + 1
      bnlocal(isame,1) = 'backbone'
      bnlocal(isame,2) = 'backbone'
      bnlocal(isame,3) = 'G'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

c     HACK FIX THIS!!!
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS_OT)
      bclocal(0) = 90.0d0
      bclocal(1) = (80.0d0)*(0.5d0)*CON_KCALTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     backbone - backbone - A [jcp2001 and pc2007]
      isame = isame + 1
      bnlocal(isame,1) = 'backbone'
      bnlocal(isame,2) = 'backbone'
      bnlocal(isame,3) = 'backbone'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * torsions
c     * cosine on diference of dihedral with equilibrium value
c     ******************************************************************

c     --- initialize tclocal
      do ii = MINTORCOEFF,MAXTORCOEFF
         tclocal(ii) = 0.0d0
      enddo

      type = type + 1
      call twh_torstyle(GLB_SET,type,23)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 80.0d0*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     base - backbone - backbone - base [jcp2001 and pc2007]
      isame = isame + 1
      tnlocal(isame,1) = 'A'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'A'
      isame = isame + 1
      tnlocal(isame,1) = 'A'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'T'
      isame = isame + 1
      tnlocal(isame,1) = 'A'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'C'
      isame = isame + 1
      tnlocal(isame,1) = 'A'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'G'
      isame = isame + 1
      tnlocal(isame,1) = 'T'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'T'
      isame = isame + 1
      tnlocal(isame,1) = 'T'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'C'
      isame = isame + 1
      tnlocal(isame,1) = 'T'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'G'
      isame = isame + 1
      tnlocal(isame,1) = 'C'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'C'
      isame = isame + 1
      tnlocal(isame,1) = 'C'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'G'
      isame = isame + 1
      tnlocal(isame,1) = 'G'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'G'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     HACK FIX THIS!!!
      type = type + 1
      call twh_torstyle(GLB_SET,type,23)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 80.0d0*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     base - backbone - backbone - backbone
      isame = isame + 1
      tnlocal(isame,1) = 'A'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'backbone'
      isame = isame + 1
      tnlocal(isame,1) = 'T'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'backbone'
      isame = isame + 1
      tnlocal(isame,1) = 'C'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'backbone'
      isame = isame + 1
      tnlocal(isame,1) = 'G'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'backbone'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     HACK FIX THIS!!!
      type = type + 1
      call twh_torstyle(GLB_SET,type,23)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 80.0d0*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     base - backbone - backbone - base [jcp2001 and pc2007]
      isame = isame + 1
      tnlocal(isame,1) = 'backbone'
      tnlocal(isame,2) = 'backbone'
      tnlocal(isame,3) = 'backbone'
      tnlocal(isame,4) = 'backbone'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
