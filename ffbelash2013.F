#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2013-2013 Marcus G. Martin                       *
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
      subroutine twh_ffbelash2013(lfailure,ntypes)
c     ******************************************************************
c     * sets up parameters for the belash2013 forcefield               *
c     *                                                                *
c     * originally written 07-16-2013 by M.G. Martin                   *
c     * last modified 07-27-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      integer ii,type,isame
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) ffname
c     --- double precision scalars
      double precision masslocal,dvalue,avalue,bvalue,cvalue
c     --- character arrays
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
c     --- double precision arrays
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)


      ffname = 'Belash2013'
      call twh_classical_potential(GLB_SET,'Embedded Atom Method')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     sets up the mercury eam potential from
c     D.K. Belashchenko, "Application of the Embedded Atom Model to 
c     Liquid Mercury", High Temperature (2013) 51, 1, pp. 40-48.

c     --- mercury
      type = 1
      masslocal = 200.59d0
      elelocal =  'Hg'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Hg'
      cbnlocal = 'Hg'
      canlocal = 'Hg'
      ctnlocal = 'Hg'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- Pair interactions with self
      ii = twh_getnbtype(type,type)
c     --- pair potential is Belashchenko 5-8
      call twh_eam_pair_style(GLB_SET,ii,'Belashchenko 5-8')
      isame = 0
c     --- r1
      isame = isame + 1
      nbclocal(isame) = 2.55d0
c     --- r2
      isame = isame + 1
      nbclocal(isame) = 3.10d0
c     --- r3
      isame = isame + 1
      nbclocal(isame) = 3.70d0
c     --- r4
      isame = isame + 1
      nbclocal(isame) = 4.60d0
c     --- r5
      isame = isame + 1
      nbclocal(isame) = 6.20d0
c     --- r6
      isame = isame + 1
      nbclocal(isame) = 8.30d0
c     --- a00
      isame = isame + 1
      nbclocal(isame) = 0.169356d0*CON_EVTOK
c     --- a01
      isame = isame + 1
      nbclocal(isame) = -6.34432d0*CON_EVTOK
c     --- a02
      isame = isame + 1
      nbclocal(isame) = 3.8d0*CON_EVTOK
c     --- a03
      isame = isame + 1
      nbclocal(isame) = 1.96d0
c     --- a10
      isame = isame + 1
      nbclocal(isame) = -0.58730732649565d-1*CON_EVTOK
c     --- a11
      isame = isame + 1
      nbclocal(isame) = -0.8725311607224d-1*CON_EVTOK
c     --- a12
      isame = isame + 1
      nbclocal(isame) = 0.23383679109512d1*CON_EVTOK
c     --- a13
      isame = isame + 1
      nbclocal(isame) = 0.54150014961337d2*CON_EVTOK
c     --- a14
      isame = isame + 1
      nbclocal(isame) = 0.50041303903954d3*CON_EVTOK
c     --- a15
      isame = isame + 1
      nbclocal(isame) = 0.22115293102511d4*CON_EVTOK
c     --- a16
      isame = isame + 1
      nbclocal(isame) = 0.50976873712314d4*CON_EVTOK
c     --- a17
      isame = isame + 1
      nbclocal(isame) = 0.58975308220171d4*CON_EVTOK
c     --- a18
      isame = isame + 1
      nbclocal(isame) = 0.27077049182679d4*CON_EVTOK
c     --- a20
      isame = isame + 1
      nbclocal(isame) = -0.55649321526289d-1*CON_EVTOK
c     --- a21
      isame = isame + 1
      nbclocal(isame) = 0.18427524715662d-1*CON_EVTOK
c     --- a22
      isame = isame + 1
      nbclocal(isame) = 0.41591385641041d0*CON_EVTOK
c     --- a23
      isame = isame + 1
      nbclocal(isame) = 0.10555632657899d2*CON_EVTOK
c     --- a24
      isame = isame + 1
      nbclocal(isame) = 0.94856510285549d2*CON_EVTOK
c     --- a25
      isame = isame + 1
      nbclocal(isame) = 0.38745647613661d3*CON_EVTOK
c     --- a26
      isame = isame + 1
      nbclocal(isame) = 0.80080059311747d3*CON_EVTOK
c     --- a27
      isame = isame + 1
      nbclocal(isame) = 0.82197252535534d3*CON_EVTOK
c     --- a28
      isame = isame + 1
      nbclocal(isame) = 0.33467726782330d3*CON_EVTOK
c     --- a30
      isame = isame + 1
      nbclocal(isame) = -0.14076511375606d-1*CON_EVTOK
c     --- a31
      isame = isame + 1
      nbclocal(isame) = 0.44487200677395d-1*CON_EVTOK
c     --- a32
      isame = isame + 1
      nbclocal(isame) = -0.85710953624001d-2*CON_EVTOK
c     --- a33
      isame = isame + 1
      nbclocal(isame) = -0.17535348275128d0*CON_EVTOK
c     --- a34
      isame = isame + 1
      nbclocal(isame) = -0.18290639200802d1*CON_EVTOK
c     --- a35
      isame = isame + 1
      nbclocal(isame) = -0.61396850926303d1*CON_EVTOK
c     --- a36
      isame = isame + 1
      nbclocal(isame) = -0.93832108328515d1*CON_EVTOK
c     --- a37
      isame = isame + 1
      nbclocal(isame) = -0.68240837932190d1*CON_EVTOK
c     --- a38
      isame = isame + 1
      nbclocal(isame) = -0.19167854338407d1*CON_EVTOK
c     --- a40
      isame = isame + 1
      nbclocal(isame) = -0.50217746756971d-2*CON_EVTOK
c     --- a41
      isame = isame + 1
      nbclocal(isame) = -0.25080999359488d-1*CON_EVTOK
c     --- a42
      isame = isame + 1
      nbclocal(isame) = -0.57924106076762d-1*CON_EVTOK
c     --- a43
      isame = isame + 1
      nbclocal(isame) = -0.21114344856450d0*CON_EVTOK
c     --- a44
      isame = isame + 1
      nbclocal(isame) = -0.63905669655741d0*CON_EVTOK
c     --- a45
      isame = isame + 1
c     --- typo in the paper (Belashchenko per comm 7-2013
      nbclocal(isame) = -0.11066271548233d1*CON_EVTOK
c     --- nbclocal(isame) = -0.11066271548233d0*CON_EVTOK
c     --- a46
      isame = isame + 1
      nbclocal(isame) = -0.10219966911147d1*CON_EVTOK
c     --- a47
      isame = isame + 1
      nbclocal(isame) = -0.46628922422428d0*CON_EVTOK
c     --- a48
      isame = isame + 1
      nbclocal(isame) = -0.82702468028215d-1*CON_EVTOK
c     --- a50
      isame = isame + 1
      nbclocal(isame) = 0.0d0*CON_EVTOK
c     --- a51
      isame = isame + 1
      nbclocal(isame) = 0.0d0*CON_EVTOK
c     --- a52
      isame = isame + 1
      nbclocal(isame) = -0.11675312540717d0*CON_EVTOK
c     --- a53
      isame = isame + 1
      nbclocal(isame) = -0.42406122320448d0*CON_EVTOK
c     --- a54
      isame = isame + 1
      nbclocal(isame) = -0.71663215147580d0*CON_EVTOK
c     --- a55
      isame = isame + 1
      nbclocal(isame) = -0.64344074447359d0*CON_EVTOK
c     --- a56
      isame = isame + 1
      nbclocal(isame) = -0.31451262110966d0*CON_EVTOK
c     --- a57
      isame = isame + 1
      nbclocal(isame) = -0.78960959177628d-1*CON_EVTOK
c     --- a58
      isame = isame + 1
      nbclocal(isame) = - 0.79477007838712d-2*CON_EVTOK
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- density with self
      call twh_eam_dens_style(GLB_SET,type,type,'exponential')
      isame = 0
c     --- density prefactor (p1 in the paper)
      isame = isame + 1
      dvalue = (4.8019d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,dvalue)
c     --- density exponential factor B (-p2 in the paper)
      isame = isame + 1
      dvalue = -(1.3095d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,dvalue)
c     --- maximum density cutoff
      dvalue = 8.30d0
      call twh_eam_dens(GLB_SET,1,type,type,isame,dvalue)
c     --- set eam_ndens
      call twh_eam_ndens(GLB_SET,type,type,isame)

c     --- Embedding function 
      call twh_eam_embed_style(GLB_SET,type,'Belashchenko 10 rho')
c     --- 
      isame = 0
c     --- isame 1 through 10 is a,b,c in slots 1,2,3
c     --- a,b,c: 1 : NOTE b1 does not exist so 0.0 spot
      isame = isame + 1
      avalue = (-0.08798d0)*CON_EVTOK
      bvalue = (0.0d0)*CON_EVTOK
      cvalue = (0.7867d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- a,b,c: 2
      isame = isame + 1
      avalue = (-0.078461d0)*CON_EVTOK
      bvalue = (-0.173074d0)*CON_EVTOK
      cvalue = (-1.40d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- a,b,c: 3
      isame = isame + 1
      avalue = (-0.073575d0)*CON_EVTOK
      bvalue = (0.050926d0)*CON_EVTOK
      cvalue = (2.00d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- a,b,c: 4
      isame = isame + 1
      avalue = (-0.058668d0)*CON_EVTOK
      bvalue = (-0.349074d0)*CON_EVTOK
      cvalue = (-3.00d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- a,b,c: 5
      isame = isame + 1
      avalue = (-0.051551d0)*CON_EVTOK
      bvalue = (0.190926d0)*CON_EVTOK
      cvalue = (5.42d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- a,b,c: 6
      isame = isame + 1
      avalue = (-0.038358d0)*CON_EVTOK
      bvalue = (-0.567874d0)*CON_EVTOK
      cvalue = (-2.00d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- a,b,c: 7
      isame = isame + 1
      avalue = (-0.005728d0)*CON_EVTOK
      bvalue = (-0.247874d0)*CON_EVTOK
      cvalue = (-3.20d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- a,b,c: 8
      isame = isame + 1
      avalue = (-0.056069d0)*CON_EVTOK
      bvalue = (0.840126d0)*CON_EVTOK
      cvalue = (4.00d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- a,b,c: 9
      isame = isame + 1
      avalue = (-0.056512d0)*CON_EVTOK
      bvalue = (0.314680d0)*CON_EVTOK
      cvalue = (0.980d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- a,b,c: 10
      isame = isame + 1
      avalue = (2.625839d0)*CON_EVTOK
      bvalue = (2.629723d0)*CON_EVTOK
      cvalue = (0.230d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- rho0,rho1,rho2
      isame = isame + 1
      avalue = (1.0d0)
      bvalue = (0.89d0)
      cvalue = (0.81d0)
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- rho3,rho4,rho5
      isame = isame + 1
      avalue = (0.71d0)
      bvalue = (0.62d0)
      cvalue = (0.55d0)
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- rho6,rho7,rho8
      isame = isame + 1
      avalue = (0.47d0)
      bvalue = (0.30d0)
      cvalue = (1.20d0)
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- rho9,m,n
      isame = isame + 1
      avalue = (2.80d0)
      bvalue = (1.70d0)
      cvalue = (3.00d0)
      call twh_eam_embed(GLB_SET,1,type,isame,avalue)
      call twh_eam_embed(GLB_SET,2,type,isame,bvalue)
      call twh_eam_embed(GLB_SET,3,type,isame,cvalue)
c     --- set nembed
      call twh_eam_nembed(GLB_SET,type,isame)


      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
