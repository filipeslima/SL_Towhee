#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2013 Marcus G. Martin                       *
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
      program twh_ffcreate
c     ******************************************************************
c     * this is the starting place for force field generation code     *
c     *                                                                *
c     * rewritten as a wrapper on twh_ffcreate_main 08-04-2011         *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- local parameters
      logical lfailure

      lfailure = .false.

      call twh_ffcreate_main(lfailure)

      if ( lfailure ) then
         write(6,*) 'FFCREATE: fatal error terminated code early'
      endif

      end

      subroutine twh_ffcreate_main(lfailure)
c     ******************************************************************
c     * controls the generation of the forcefield files                *
c     *                                                                *
c     * originally written 06-21-2002 by M.G. Martin                   *
c     * last modified 07-19-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBSAME
#define FUNCTION_GET_NTSAME
#define FUNCTION_GET_NVSAME
#include "functions.h"
c     --- values passed to/from the subroutine
      logical lfailure
c     --- local parameters
      integer maxff
      parameter (maxff=100)
c     --- local variables
c     --- character string scalars
      character*(FFNAMELEN) forcefield
      character*(MAXDIRLENGTH) ffname
c     --- character string arrays
      character*(FFNAMELEN) fflist
      dimension fflist(maxff)
c     --- logical scalars
      logical lall
c     --- logical arrays
      logical ffoutput
      dimension ffoutput(maxff)
c     --- integer scalars
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
      integer naatypes,nbitypes,noftypes,nhbondtypes
      integer iffnum,nff,iff,istart,istop,type
      integer maxtypes,maxbondtypes,maxangletypes,maxdihedtypes
      integer maximprotypes,maxaatypes,maxbitypes,maxoftypes
      integer maxsamevib,maxsameben,maxsametor,maxhbondtypes
      integer ffcount

c     --- initialize max variables
      maxtypes = 0
      maxbondtypes = 0
      maxangletypes = 0
      maxdihedtypes = 0
      maximprotypes = 0
      maxaatypes = 0
      maxbitypes = 0
      maxoftypes = 0
      maxhbondtypes = 0
      ffcount = 0
c     --- initialize maxsame variables
      maxsamevib = 0
      maxsameben = 0
      maxsametor = 0

c     --- here is the complete list of force fields along with a logical
c     --- that determines whether this is part of those generated with
c     --- the "all" command
c     --- those not automatically generated are ones that require file
c     --- translation

c     --- initialize nff
      nff = 0
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Ackl2004'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Alavi2005'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Amber96'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'Aqvist'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'Belash2006'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'Belash2013'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'C19eef1'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'C19sasa'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'C27rigid'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'CatlowFaux'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'Charmm19'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'Charmm22'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'Charmm22fe'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'Charmm27'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'Charmm27x'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'CharmmFile'
      ffoutput(nff) = .false.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) =  'ClayFF'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'COMPASSv1'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Coon1987'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Cui1998'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Cui2002'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'DACNIS-UA'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'devel'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Dick1994'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Ding1986'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'DREIDING'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Dubb2004'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'DWS2001'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Elli2002'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'EPM'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Fris2003'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Fris2008'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Gala1994'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Gordon'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Gromos43A1'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'HardSphere'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Jaramillo'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'KFvBvS'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Last1993'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'LGM'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'LJium'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Mend2003'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'mgmstereo'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'MCY1976'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'MM2'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'MMFF94'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Morrow2002'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'NERDv1'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'NERDv2'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'NERDv3'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'OPLS-1996'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'OPLS-2001'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'OPLS-aa'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'OPLS-ua'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Pana1989'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'PMF'
      ffoutput(nff) = .false.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Potter1997'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'QMFF-VIII'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Richar1995'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Setfl'
      ffoutput(nff) = .false.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Shah2004'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Shukla1987'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'SKS'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Smith1994'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'SMMKmain'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'SMMKnaip'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'SPC-E'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'SquareWell'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Stil1985'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Sum2003'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Tele1987'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'TIP3P'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'TIP4P'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'TIP5P'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'TraPPE-EH'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'TraPPE-UA'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'TraPPE-UAf'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'UFF'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Unlu2004'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Vega1992'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Vink2001'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Walt2001'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Weiner1984'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Weiner1986'
      ffoutput(nff) = .true.
c     --- add a forcefield to the list
      nff = nff + 1
      fflist(nff) = 'Wielop1985'
      ffoutput(nff) = .true.

c     --- set the forcefield
      write(6,*) 'Choose the name of the forcefield to generate'
      write(6,*) 'a Towhee forcefield file.'
      write(6,*) 'All'
      do iff = 1,nff
         write(6,*) fflist(iff)
      enddo
      read(5,*) forcefield

      if ( forcefield .eq. 'All' ) then
         lall = .true.
         istart = 1
         istop = nff
      else
         lall = .false.
         istart = 1
         istop = 1
      endif

      do iff = istart,istop
         call twh_initialize(lfailure)
         if ( lfailure ) then
            write(6,*) 'FFCREATE_MAIN: fatal error in initialze'
            return
         endif

         if ( lall ) then
c           --- kick out if this is a file translator
            if ( .not. ffoutput(iff) ) goto 200
c           --- set the forcefield
            forcefield = fflist(iff)
         endif

c        --- initialize the counters
         ntypes = 0
         nbondtypes = 0
         nangletypes = 0
         ndihedtypes = 0
         nimprotypes = 0
         naatypes = 0
         nbitypes = 0
         noftypes = 0
         nhbondtypes = 0

         if ( forcefield .eq. 'Ackl2004' ) then
            call twh_ffackl2004(lfailure,ntypes)
         elseif ( forcefield .eq. 'Alavi2005' ) then
            call twh_ffalavi2005(lfailure,ntypes,nbondtypes,nangletypes
     &           ,nbitypes)
         elseif ( forcefield .eq. 'Amber96' ) then
            call twh_ffamber96(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'Aqvist' ) then
            call twh_ffaqvist(lfailure,ntypes)
         elseif ( forcefield .eq. 'Belash2006' ) then
            call twh_ffbelash2006(lfailure,ntypes)
         elseif ( forcefield .eq. 'Belash2013' ) then
            call twh_ffbelash2013(lfailure,ntypes)
         elseif ( forcefield .eq. 'CatlowFaux' ) then
            call twh_ffcatlowfaux(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes)
         elseif ( forcefield .eq. 'Charmm19' ) then
            call twh_ffcharmm19(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'C19eef1' ) then
            call twh_ffc19eef1(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'C19sasa' ) then
            call twh_ffc19sasa(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'C27rigid' ) then
            call twh_ffc27rigid(lfailure,ntypes,nbondtypes,nangletypes
     &           ,nbitypes)
         elseif ( forcefield .eq. 'Charmm22' ) then
            call twh_ffcharmm22(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'Charmm22fe' ) then
            call twh_ffcharmm22fe(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes)
         elseif ( forcefield .eq. 'Charmm27' ) then
            call twh_ffcharmm27(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'Charmm27x' ) then
            call twh_ffcharmm27x(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'CharmmFile') then
            call twh_ffreadcharmmfile(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes,nimprotypes,ffname)
         elseif ( forcefield .eq. 'ClayFF' ) then
            call twh_ffclayff(lfailure,ntypes,nbondtypes,nangletypes)
         elseif ( forcefield .eq. 'COMPASSv1' ) then
            call twh_ffcompassv1(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes,naatypes,nbitypes)
         elseif ( forcefield .eq. 'Coon1987' ) then
            call twh_ffcoon1987(lfailure,ntypes,nbondtypes,nbitypes)
         elseif ( forcefield .eq. 'Cui1998' ) then
            call twh_ffcui1998(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nbitypes)
         elseif ( forcefield .eq. 'Cui2002' ) then
            call twh_ffcui2002(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'DACNIS-UA' ) then
            call twh_ffdacnisua(lfailure,ntypes,nbondtypes)
         elseif ( forcefield .eq. 'devel' ) then
            call twh_ffdevel(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'Dick1994' ) then
            call twh_ffdick1994(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'Ding1986' ) then
            call twh_ffding1986(lfailure,ntypes)
         elseif ( forcefield .eq. 'DREIDING' ) then
            call twh_ffdreiding(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'Dubb2004' ) then
            call twh_ffdubb2004(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nbitypes)
         elseif ( forcefield .eq. 'DWS2001' ) then
            call twh_ffdws2001(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'Elli2002' ) then
            call twh_ffelli2002(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'EPM' ) then
            call twh_ffepm(lfailure,ntypes,nbondtypes,nangletypes
     &           ,nbitypes)
         elseif ( forcefield .eq. 'Fris2003' ) then
            call twh_fffris2003(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'Fris2008' ) then
            call twh_fffris2008(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'Gala1994' ) then
            call twh_ffgala1994(lfailure,ntypes,nbondtypes,nbitypes)
         elseif ( forcefield .eq. 'Gordon' ) then
            call twh_ffgordon(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'Gromos43A1' ) then
            call twh_ffgromos43a1(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'HardSphere' ) then
            call twh_ffhardsphere(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes)
         elseif ( forcefield .eq. 'Jaramillo' ) then
            call twh_ffjaramillo(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'KFvBvS' ) then
            call twh_ffkfvbvs(lfailure,ntypes)
         elseif ( forcefield .eq. 'Last1993' ) then
            call twh_fflast1993(lfailure,ntypes)
         elseif ( forcefield .eq. 'LGM' ) then
            call twh_fflgm(lfailure,ntypes)
         elseif ( forcefield .eq. 'LJium' ) then
            call twh_ffljium(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'MCY1976' ) then
            call twh_ffmcy1976(lfailure,ntypes,nbondtypes,nangletypes
     &           ,nbitypes)
         elseif ( forcefield .eq. 'Mend2003' ) then
            call twh_ffmend2003(lfailure,ntypes)
         elseif ( forcefield .eq. 'mgmstereo' ) then
            call twh_ffmgmstereo(lfailure,nimprotypes)
         elseif ( forcefield .eq. 'MM2' ) then
            call twh_ffmm2(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'MMFF94' ) then
            call twh_ffmmff94(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes,nbitypes)
         elseif ( forcefield .eq. 'Morrow2002' ) then
            call twh_ffmorrow2002(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes,nimprotypes,nbitypes)
         elseif ( forcefield .eq. 'NERDv1' ) then
            call twh_ffnerdv1(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'NERDv2' ) then
            call twh_ffnerdv2(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'NERDv3' ) then
            call twh_ffnerdv3(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'OPLS-aa' ) then
            call twh_ffoplsaa(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'OPLS-ua' ) then
            call twh_ffoplsua(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'OPLS-1996' ) then
            call twh_ffopls1996(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'OPLS-2001' ) then
            call twh_ffopls2001(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'Pana1989' ) then
            call twh_ffpana1989(lfailure,ntypes)
         elseif ( forcefield .eq. 'PMF' ) then
            call twh_ffpmf(lfailure,ntypes)
         elseif ( forcefield .eq. 'Potter1997' ) then
            call twh_ffpotter1997(lfailure,ntypes,nbondtypes
     &           ,nangletypes,nbitypes)
         elseif ( forcefield .eq. 'QMFF-VIII' ) then
            call twh_ffqmff_viii(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,naatypes,nbitypes)
         elseif ( forcefield .eq. 'Setfl' ) then
            call twh_ffreadsetflfile(lfailure,ntypes,forcefield)
         elseif ( forcefield .eq. 'Richar1995' ) then
            call twh_ffrichar1995(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes,nbitypes)
         elseif ( forcefield .eq. 'Shah2004' ) then
            call twh_ffshah2004(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nbitypes)
         elseif ( forcefield .eq. 'Shukla1987' ) then
            call twh_ffshukla1987(lfailure,ntypes)
         elseif ( forcefield .eq. 'SKS' ) then
            call twh_ffsks(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'Smith1994' ) then
            call twh_ffsmith1994(lfailure,ntypes)
         elseif ( forcefield .eq. 'SMMKmain' ) then
            call twh_ffsmmkmain(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,noftypes)
         elseif ( forcefield .eq. 'SMMKnaip' ) then
            call twh_ffsmmknaip(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'SPC-E' ) then
            call twh_ffspce(lfailure,ntypes,nbondtypes,nangletypes)
         elseif ( forcefield .eq. 'SquareWell' ) then
            call twh_ffsquarewell(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes)
         elseif ( forcefield .eq. 'Stil1985' ) then
            call twh_ffstil1985(lfailure,ntypes)
         elseif ( forcefield .eq. 'Sum2003' ) then
            call twh_ffsum2003(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'Tele1987' ) then
            call twh_fftele1987(lfailure,ntypes,nbondtypes,nangletypes
     &           ,nbitypes)
         elseif ( forcefield .eq. 'TIP3P' ) then
            call twh_fftip3p(lfailure,ntypes,nbondtypes,nangletypes
     &           ,nbitypes)
         elseif ( forcefield .eq. 'TIP4P' ) then
            call twh_fftip4p(lfailure,ntypes,nbondtypes,nangletypes
     &           ,nimprotypes,noftypes)
         elseif ( forcefield .eq. 'TIP5P' ) then
            call twh_fftip5p(lfailure,ntypes,nbondtypes,nangletypes)
         elseif ( forcefield .eq. 'TraPPE-EH' ) then
            call twh_fftrappeeh(lfailure,ntypes,nbondtypes,nangletypes
     &           ,nbitypes)
         elseif ( forcefield .eq. 'TraPPE-UA' ) then
            call twh_fftrappeua(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nbitypes,noftypes)
         elseif ( forcefield .eq. 'TraPPE-UAf' ) then
            call twh_fftrappeuaf(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nbitypes)
         elseif ( forcefield .eq. 'UFF' ) then
            call twh_ffuff(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nimprotypes)
         elseif ( forcefield .eq. 'Unlu2004' ) then
            call twh_ffunlu2004(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'Vega1992' ) then
            call twh_ffvega1992(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes)
         elseif ( forcefield .eq. 'Vink2001' ) then
            call twh_ffvink2001(lfailure,ntypes)
         elseif ( forcefield .eq. 'Walt2001' ) then
            call twh_ffwalt2001(lfailure,ntypes,nbondtypes,nangletypes
     &           ,ndihedtypes,nbitypes)
         elseif ( forcefield .eq. 'Weiner1984' ) then
            call twh_ffweiner1984(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes,nimprotypes,nhbondtypes)
         elseif ( forcefield .eq. 'Weiner1986' ) then
            call twh_ffweiner1986(lfailure,ntypes,nbondtypes
     &           ,nangletypes,ndihedtypes,nimprotypes,nhbondtypes)
         elseif ( forcefield .eq. 'Wielop1985' ) then
            call twh_ffwielop1985(lfailure,ntypes,nbondtypes
     &           ,nangletypes,nbitypes)
         else   
            write(6,*) 'unknown force field'
            write(6,*) forcefield
            lfailure = .true.
            return
         endif
         if ( lfailure ) return

c        --- set the output
         iffnum = 1
         ffname = 'towhee_ff_'
         ffname(11:20) = forcefield
         call twh_ff_filename(GLB_SET,iffnum,ffname)

         call twh_rwforcefield(lfailure,WRITE_FLAG,ntypes,nbondtypes
     &        ,nangletypes,ndihedtypes,nimprotypes,naatypes,nbitypes
     &        ,noftypes,nhbondtypes,iffnum)
         if ( lfailure ) then
            write(6,*) 'FFCREATE_MAIN: fatal error in rwforcefield'
            return
         endif
         ffcount = ffcount + 1

         write(6,*) 'A force field file for the ',forcefield
     &        ,' forcefield ','was written to ',ffname

c        --- update maximums
         if ( ntypes .gt. maxtypes ) maxtypes = ntypes
         if ( nbondtypes .gt. maxbondtypes ) maxbondtypes = nbondtypes
         if ( nangletypes .gt. maxangletypes ) 
     &        maxangletypes = nangletypes
         if ( ndihedtypes .gt. maxdihedtypes ) 
     &        maxdihedtypes = ndihedtypes
         if ( nimprotypes .gt. maximprotypes ) 
     &        maximprotypes = nimprotypes
         if ( naatypes .gt. maxaatypes ) maxaatypes = naatypes
         if ( nbitypes .gt. maxbitypes ) maxbitypes = nbitypes
         if ( noftypes .gt. maxoftypes ) maxoftypes = noftypes
         if ( nhbondtypes .gt. maxhbondtypes )
     &        maxhbondtypes = nhbondtypes
c        --- update maxsames
         do type = 1,nbondtypes
            if ( twh_get_nvsame(type) .gt. maxsamevib ) then
               maxsamevib = twh_get_nvsame(type)
            endif
         enddo
         do type = 1,nangletypes
            if ( twh_get_nbsame(type) .gt. maxsameben ) then
               maxsameben = twh_get_nbsame(type)
            endif
         enddo
         do type = 1,ndihedtypes
            if ( twh_get_ntsame(type) .gt. maxsametor ) then
               maxsametor = twh_get_ntsame(type)
            endif
         enddo
 200     continue
      enddo

      if ( istop .gt. 1 ) then
         write(6,*) 'Number of forcefield files created:',ffcount
         write(6,*) 'Maximum Number of Each Type for preproc.h'
         write(6,*) 'types, NNTYPE:',maxtypes,NNTYPE
         write(6,*) 'bonds, TVIBMAX:',maxbondtypes,TVIBMAX
         write(6,*) 'angles, TBENMAX:',maxangletypes,TBENMAX
         write(6,*) 'dihed, TTORMAX:',maxdihedtypes,TTORMAX
         write(6,*) 'improp, TIMPMAX:',maximprotypes,TIMPMAX
         write(6,*) 'aa, TAAMAX:',maxaatypes,TAAMAX
         write(6,*) 'bond incr TBIMAX:',maxbitypes,TBIMAX
         write(6,*) 'one-five TOFMAX:',maxoftypes,TOFMAX
         write(6,*) 'hbond, THBONDMAX:',maxhbondtypes,THBONDMAX
         write(6,*) 'vsame MAXSAMEV',maxsamevib,MAXSAMEV
         write(6,*) 'bsame MAXSAMEB',maxsameben,MAXSAMEB
         write(6,*) 'tsame MAXSAMET',maxsametor,MAXSAMET
      endif

      return
      end

      subroutine twh_simple_distance(dist,apos,bpos)
c     ******************************************************************
c     * computes a simple distance between two points without any      *
c     * concern for minimum image conventions                          *
c     *                                                                *
c     * originally written 04-18-2006 by M.G. Martin                   *
c     * last modified 04-18-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision dist,apos,bpos
      dimension apos(3),bpos(3)
c     --- local variable
      integer icord

      dist = 0.0d0
      do icord = 1,3
         dist = dist + (apos(icord)-bpos(icord))**2
      enddo
      dist = dsqrt(dist)
      return
      end

      subroutine twh_simple_angle(angle,xcord,ycord,zcord)
c     *****************************************************************
c     * computes an angle from coordinates without any minimum image  *
c     *                                                               *
c     * originally written 04-20-2006 by M.G. Martin                  *
c     * last modified 04-20-2006 by M.G. Martin                       *
c     *****************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_DOTPRODUCT
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision angle,xcord,ycord,zcord
      dimension xcord(3),ycord(3),zcord(3)
c     --- local variables
      double precision xone,yone,zone,xtwo,ytwo,ztwo,distone,disttwo
      double precision angcos

      distone = 0.0d0
      disttwo = 0.0d0

      xone = xcord(1)-xcord(2)
      yone = ycord(1)-ycord(2)
      zone = zcord(1)-zcord(2)
      xtwo = xcord(3)-xcord(2)
      ytwo = ycord(3)-ycord(2)
      ztwo = zcord(3)-zcord(2)

      angcos = twh_dotproduct(xone,yone,zone,xtwo,ytwo,ztwo,distone
     &     ,disttwo)
      angle = twh_arccos(angcos)
      return
      end

      subroutine twh_simple_dihedral(dihedral,xcord,ycord,zcord)
c     ******************************************************************
c     * computes a dihedral from coordinates without any minimum image *
c     *                                                                *
c     * originally written 04-25-2006 by M.G. Martin                   *
c     * last modified 05-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_CROSSPRODUCT
#define FUNCTION_DOTPRODUCT
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision dihedral,xcord,ycord,zcord
      dimension xcord(4),ycord(4),zcord(4)
c     --- local variables
      integer ivec
      double precision xvec,yvec,zvec
      dimension xvec(3),yvec(3),zvec(3)
      double precision xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2,phic

      do ivec = 1,3
         xvec(ivec) = xcord(ivec+1) - xcord(ivec)
         yvec(ivec) = ycord(ivec+1) - ycord(ivec)
         zvec(ivec) = zcord(ivec+1) - zcord(ivec)
      enddo
c     --- calculate cross products 1-2 x 2-3 
      call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(2),yvec(2)
     &     ,zvec(2),xaa1,yaa1,zaa1)
c     --- calculate cross products 3-4 x 2-3
      call twh_crossproduct(xvec(3),yvec(3),zvec(3),xvec(2),yvec(2)
     &     ,zvec(2),xaa2,yaa2,zaa2)
c     --- calculate dot product of cross products
      phic = -twh_dotproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2)
      dihedral = twh_arccos(phic)

      return
      end
