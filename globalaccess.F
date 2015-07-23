#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2008 Marcus G. Martin                       *
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

c     ******************************************************************
c     * this file contains all of the functions and special subroutines*
c     * for accessing global data along with subroutines that bundle   *
c     * the adjustment of multiple global data structures              *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************

      subroutine twh_create_torsion(imolty,iunit,itor2,itor3,itor4
     &     ,itype)
c     ******************************************************************
c     * creates a new torsion by incrementing the number of torsions,  *
c     * assigning the values of ijtor2, ijtor3, ijtor4; and setting the*
c     * new torsion type                                               *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer imolty,iunit,itor2,itor3,itor4,itype
c     --- local variables
      integer itemp

      call twh_intor(GLB_INCR,imolty,iunit,1)
      call twh_intor(GLB_GET,imolty,iunit,itemp)
      if ( itemp .gt. MAXTOR ) then
         write(6,*) 'CREATE_TORSION: new torsion exceeds MAXTOR'
         stop
      endif
      call twh_ijtor2(GLB_SET,imolty,iunit,itemp,itor2)
      call twh_ijtor3(GLB_SET,imolty,iunit,itemp,itor3)
      call twh_ijtor4(GLB_SET,imolty,iunit,itemp,itor4)
      call twh_ittor(GLB_SET,imolty,iunit,itemp,itype)

      return
      end

      subroutine twh_set_aacoeff(itype,coeffs)
c     ******************************************************************
c     * assigns the aacoeff variables for this type                    *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffs
      dimension coeffs(MINAACOEFF:MAXAACOEFF)
c     --- local variables
      integer index

      do index = MINAACOEFF,MAXAACOEFF
         call twh_aacoeff(GLB_SET,itype,index,coeffs(index))
c        --- reset the coeff back to zero for the next assignment
         coeffs(index) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_aanames(itype,nsame,namearray)
c     ******************************************************************
c     * sets the number of quartets with the same parameters and also  *
c     * set the names of those quartets                                *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMEAA,4)
c     --- local variables
      integer isame,index

      call twh_naasame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do index = 1,4
            call twh_aanames(GLB_SET,itype,isame,index
     &           ,namearray(isame,index))
c           --- reset for next assignment
            namearray(isame,index) = '-'
         enddo
      enddo
      nsame = 0

      return
      end

      subroutine twh_set_bencoeff(itype,coeffs)
c     ******************************************************************
c     * assigns the bencoeff variables for this type                   *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffs
      dimension coeffs(MINBENDCOEFF:MAXBENDCOEFF)
c     --- local variables
      integer index

      do index = MINBENDCOEFF,MAXBENDCOEFF
         call twh_bencoeff(GLB_SET,itype,index,coeffs(index))
c        --- reset the coeff back to zero for the next assignment
         coeffs(index) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_bendnames(itype,nsame,namearray)
c     ******************************************************************
c     * assigns all of the bendnames and sets nbsame for this type of  *
c     * bending angle                                                  *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMEB,MAXBNAME)
c     --- local variables
      integer isame,iatom

      call twh_nbsame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do iatom = 1,MAXBNAME
            call twh_bendnames(GLB_SET,itype,isame,iatom
     &           ,namearray(isame,iatom))
c           --- now reset the namearray for the next allocation
            namearray(isame,iatom) = '-'
         enddo
      enddo

      return
      end

      subroutine twh_set_classicalpotnum(cplocal)
c     ******************************************************************
c     * translates the classical_potential string definition into a    *
c     * defined integer to speed the inner loop in vtwobody            *
c     *                                                                *
c     * split from A.A. Chen implementation in readclassical 11-10-2005*
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) cplocal
c     --- local variables
      integer classical_pot_num

c     --- these definitions must match those in preproc.h
c     --- note all tabular potentials are numbered > 100
      if (cplocal .eq. '12-6 plus 12-10 H-bond') then
         classical_pot_num = CP_12_6_PLUS_12_10_H_BOND
      elseif (cplocal .eq. '12-6 plus solvation') then
         classical_pot_num = CP_12_6_PLUS_SOLVATION
      elseif (cplocal .eq. '12-9-6') then
         classical_pot_num = CP_12_9_6
      elseif ( cplocal .eq. 'Hard 2580 Multistep' ) then
         classical_pot_num = CP_HARD_2580_MULTISTEP
      elseif (cplocal .eq. '9-6') then
         classical_pot_num = CP_9_6
      elseif ( cplocal .eq. 'Buffered 14-7' ) then
         classical_pot_num = CP_BUFFERED_14_7
      elseif ( cplocal .eq. 'Double Exponential' ) then
         classical_pot_num = CP_DOUBLE_EXPONENTIAL
      elseif ( cplocal .eq. 'Drukker H-bond' ) then
         classical_pot_num = CP_DRUKKER_H_BOND
      elseif (cplocal .eq. 'EAM pair only' ) then
         classical_pot_num = CP_EAM_PAIR_ONLY
      elseif (cplocal .eq. 'Embedded Atom Method') then
         classical_pot_num = CP_EMBEDDED_ATOM_METHOD
      elseif (cplocal .eq. 'Exponential-12-6') then
         classical_pot_num = CP_EXPONENTIAL_12_6
      elseif (cplocal .eq. 'Exponential-6') then
         classical_pot_num = CP_EXPONENTIAL_6
      elseif (cplocal .eq. 'Gordon n-6') then
         classical_pot_num = CP_GORDON_N_6
      elseif (cplocal .eq. 'Hard Sphere') then
         classical_pot_num = CP_HARD_SPHERE
      elseif (cplocal .eq. 'Lennard-Jones') then
         classical_pot_num = CP_LENNARD_JONES
      elseif (cplocal .eq. 'Multiwell') then
         classical_pot_num = CP_MULTIWELL
      elseif ( cplocal .eq. 'Repulsive 2580 Multistep' 
     &        ) then
         classical_pot_num = CP_REPULSIVE_2580_MULTISTEP
      elseif (cplocal .eq. 'Repulsive Multiwell') then
         classical_pot_num = CP_REPULSIVE_MULTIWELL
      elseif (cplocal .eq. 'Repulsive Sphere') then
         classical_pot_num = CP_REPULSIVE_SPHERE
      elseif (cplocal .eq. 'Repulsive Well') then
         classical_pot_num = CP_REPULSIVE_WELL
      elseif (cplocal .eq. 'Scaled Lennard-Jones') then
         classical_pot_num = CP_SCALED_LENNARD_JONES
      elseif (cplocal .eq. 'Square Well') then
         classical_pot_num = CP_SQUARE_WELL
      elseif (cplocal .eq. 'Stillinger-Weber') then
         classical_pot_num = CP_STILLINGER_WEBER
      elseif ( cplocal .eq. 'SW pair only' ) then
         classical_pot_num = CP_SW_PAIR_ONLY
      elseif (cplocal .eq. 'Tabulated Pair') then
         classical_pot_num = CP_TABULATED_PAIR
      elseif (cplocal .eq. 'UFF 12-6') then
         classical_pot_num = CP_UFF_12_6
      elseif (cplocal .eq. 'Weeks-Chandler-Anderson') then
         classical_pot_num = CP_WEEKS_CHANDLER_ANDERSON
      elseif ( cplocal .eq. 'init-error' .or. cplocal .eq. 'wild' ) then
c        --- this value is set as the initial classical_potential to 
c        --- detect errors
         classical_pot_num = CP_INIT_ERROR
      else 
         write(6,*) 'SETCLASSICALPOTNUM: invalid classical_potential'
         write(6,'(a)') cplocal
         stop
      endif

      call twh_cp_number(GLB_SET,classical_pot_num)

      return
      end

      subroutine twh_set_eamdensnum(clocal,ipair,jpair)
c     ******************************************************************
c     * translates the eam_dens_style into a defined integer to speed  *
c     * the compares                                                   *
c     *                                                                *
c     * originally written 07-20-2013 by M.G. Martin                   *
c     * last modified 07-21-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) clocal
      integer ipair,jpair
c     --- local variables
      integer eam_dens_num

c     --- these definitions must match those in preproc.h
c     --- note all tabular potentials are numbered > 100
      if ( clocal .eq. 'Ackland cubic sum' ) then
         eam_dens_num = EAM_DENS_ACKLAND
      elseif ( clocal .eq. 'exponential' ) then
         eam_dens_num = EAM_DENS_EXPON
      elseif ( clocal .eq. 'table' ) then
         eam_dens_num = EAM_DENS_TABLE
      elseif ( clocal .eq. 'none' ) then
         eam_dens_num = EAM_DENS_NONE
      else 
         write(6,*) 'SET_EAMDENSNUM: invalid eamdensnum'
         write(6,'(a)') clocal
         stop
      endif

      call twh_eam_dens_num(GLB_SET,ipair,jpair,eam_dens_num)

      return
      end

      subroutine twh_set_eamembednum(clocal,ipair)
c     ******************************************************************
c     * translates the eam_embed_style into a defined integer to speed *
c     * the compares                                                   *
c     *                                                                *
c     * originally written 07-21-2013 by M.G. Martin                   *
c     * last modified 07-21-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) clocal
      integer ipair
c     --- local variables
      integer itest

c     --- these definitions must match those in preproc.h
c     --- note all tabular potentials are numbered > 100
      if ( clocal .eq. 'Belashchenko' ) then
         itest = EAM_EMBED_BELASH
      elseif ( clocal .eq. 'Belashchenko 10 rho' ) then
         itest = EAM_EMBED_BELASHTRHO
      elseif ( clocal .eq. 'logarithmic' ) then
         itest = EAM_EMBED_LOGARITHMIC
      elseif ( clocal .eq. 'power 0.5 and 2' ) then
         itest = EAM_EMBED_POWHT
      elseif ( clocal .eq. 'power 0.5, 2, and 4' ) then
         itest = EAM_EMBED_POWHTF
      elseif ( clocal .eq. 'squareroot' ) then
         itest = EAM_EMBED_SR
      elseif ( clocal .eq. 'table' ) then
         itest = EAM_EMBED_TABLE
      else 
         write(6,*) 'SET_EAMEMBEDNUM: invalid eamembedstyle'
         write(6,'(a)') clocal
         stop
      endif

      call twh_eam_embed_num(GLB_SET,ipair,itest)

      return
      end

      subroutine twh_set_eampairnum(clocal,ipair)
c     ******************************************************************
c     * translates the eam_pair_style into a defined integer to speed  *
c     * the compares                                                   *
c     *                                                                *
c     * originally written 07-20-2013 by M.G. Martin                   *
c     * last modified 07-21-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) clocal
      integer ipair
c     --- local variables
      integer eam_pair_num

c     --- these definitions must match those in preproc.h
c     --- note all tabular potentials are numbered > 100
      if ( clocal .eq. 'Ackland 3-part' ) then
         eam_pair_num = EAM_PAIR_ACKLANDTP
      elseif ( clocal .eq. 'Ackland Power' ) then
         eam_pair_num = EAM_PAIR_ACKLANDPOW
      elseif ( clocal .eq. 'exponential' ) then
         eam_pair_num = EAM_PAIR_EXPON
      elseif ( clocal .eq. 'morse' ) then
         eam_pair_num = EAM_PAIR_MORSE
      elseif ( clocal .eq. 'table' ) then
         eam_pair_num = EAM_PAIR_TABLE
      elseif ( clocal .eq. 'Belashchenko 5-8' ) then
         eam_pair_num = EAM_PAIR_BELASHFE
      elseif ( clocal .eq. 'none' ) then
         eam_pair_num = EAM_PAIR_NONE
      else 
         write(6,*) 'SET_EAMPAIRNUM: invalid eampairnum'
         write(6,'(a)') clocal
         stop
      endif

      call twh_eam_pair_number(GLB_SET,ipair,eam_pair_num)

      return
      end

      subroutine twh_set_hbondcoeff(itype,coeffarray)
c     ******************************************************************
c     * sets up all of the coefficients for a single special hbond type*
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffarray
      dimension coeffarray(MINHBONDCOEFF:MAXHBONDCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINHBONDCOEFF,MAXHBONDCOEFF
         call twh_hbondcoeff(GLB_SET,itype,icoeff,coeffarray(icoeff))
c        --- now reset the coeffarray for the next allocation
         coeffarray(icoeff) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_hbondnames(itype,namearray)
c     ******************************************************************
c     * sets up all of the atom names for a single special hbond type  *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      character*(FFNAMELEN) namearray
      dimension namearray(MAXHBONDNAME)
c     --- local variables
      integer icoeff

      do icoeff = 1,2
         call twh_hbondnames(GLB_SET,itype,icoeff,namearray(icoeff))
c        --- now reset the namearray for the next allocation
         namearray(icoeff) = '-'
      enddo

      return
      end

      subroutine twh_set_impcoeff(itype,coeffarray)
c     ******************************************************************
c     * assigns all of the improper torsion coefficients and then      *
c     * initializes the coeffarray for the next batch                  *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffarray
      dimension coeffarray(MINIMPCOEFF:MAXIMPCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINIMPCOEFF,MAXIMPCOEFF
         call twh_impcoeff(GLB_SET,itype,icoeff,coeffarray(icoeff))
c        --- now reset the coeffarray for the next allocation
         coeffarray(icoeff) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_impnames(itype,nsame,namearray)
c     ******************************************************************
c     * assigns all of the impnames and sets nimpsame for this type of *
c     * improper torsion angle                                         *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMEIMP,MAXIMPNAME)
c     --- local variables
      integer isame,iatom

      call twh_nimpsame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do iatom = 1,MAXIMPNAME
            call twh_impnames(GLB_SET,itype,isame,iatom
     &           ,namearray(isame,iatom))
c           --- now reset the namearray for the next allocation
            namearray(isame,iatom) = '-'
         enddo
      enddo

      return
      end

      subroutine twh_set_isolvtype(lfailure,solvation_style
     &     ,solvation_type)
c     ******************************************************************
c     * sets the integer code for isolvtype based on the               *
c     * solvation_style and the solvation_type                         *
c     *                                                                *
c     * originally written 01-12-2007 by M.G. Martin                   *
c     * last modified 01-12-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*(*) solvation_style,solvation_type
c     --- local variables
      integer isolv_type

      if ( solvation_style .eq. 'none' ) then
         isolv_type = SOLV_NONE
      else
         if ( solvation_type .eq. 'tramonto' ) then
            isolv_type = SOLV_TRAMONTO
         elseif ( solvation_type .eq. 'SASA' ) then
            isolv_type = SOLV_SASA
         elseif ( solvation_type .eq. 'EEF1' ) then
            isolv_type = SOLV_EEF1
         else
            write(6,*) 'SET_ISOLVTYPE: unknown solvation_type'
            lfailure = .true.
            return
         endif
      endif

      call twh_isolvtype(GLB_SET,isolv_type)

      return
      end

      subroutine twh_set_names(itype,nbname,bondname,anglename
     &     ,torsionname)
c     ******************************************************************
c     * assigns all of the names for each nonbonded atom type for use  *
c     * in the automatic assignment of parameters                      *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-17-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      character*(*) nbname,bondname,anglename,torsionname

c     --- safety check type
      if ( itype .gt. NNTYPE ) then
         write(6,*) 'insufficient memory for this nonbond assignment'
         write(6,*) 'requested element:',itype
         write(6,*) 'NNTYPE:',NNTYPE
         stop
      endif

      call twh_nbname(GLB_SET,itype,nbname)
      call twh_cbname(GLB_SET,itype,bondname)
      call twh_caname(GLB_SET,itype,anglename)
      call twh_ctname(GLB_SET,itype,torsionname)
c     --- reset all of the names to error values
      nbname = 'init-error'
      bondname = 'init-error'
      anglename = 'init-error'
      torsionname = 'init-error'

      return
      end

      subroutine twh_set_ofnames(itype,namearray)
c     ******************************************************************
c     * assigns all of the one-five names for a single type            *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      character*(FFNAMELEN) namearray
      dimension namearray(MAXOFNAME)
c     --- local variables
      integer icoeff

      do icoeff = 1,MAXOFNAME
         call twh_ofnames(GLB_SET,itype,icoeff,namearray(icoeff))
c        --- reset all of the names to error values
         namearray(icoeff) = '-'
      enddo

      return
      end

      subroutine twh_set_qbinames(itype,namearray)
c     ******************************************************************
c     * assigns all of the qbinames for this type                      *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      character*(FFNAMELEN) namearray
      dimension namearray(2)
c     --- local variables
      integer iatom

      do iatom = 1,2
         call twh_qbinames(GLB_SET,itype,iatom,namearray(iatom))
c        --- now reset the namearray for the next allocation
         namearray(iatom) = '-'
      enddo

      return
      end

      subroutine twh_set_nbcoeff(itype,coeffarray)
c     ******************************************************************
c     * assigns all of the nonbonded coefficients for this cross type  *
c     *                                                                *
c     * originally written 04-11-2006 by M.G. Martin                   *
c     * last modified 04-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffarray
      dimension coeffarray(MINNBCOEFF:MAXNBCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINNBCOEFF,MAXNBCOEFF
         call twh_nbcoeff(GLB_SET,icoeff,itype,coeffarray(icoeff))
c        --- reset for the next assignments
         coeffarray(icoeff) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_torcoeff(itype,coeffarray)
c     ******************************************************************
c     * assigns all of the torsion coefficients and then initializes   *
c     * the coeffarray for the next batch                              *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffarray
      dimension coeffarray(MINTORCOEFF:MAXTORCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINTORCOEFF,MAXTORCOEFF
         call twh_torcoeff(GLB_SET,itype,icoeff,coeffarray(icoeff))
c        --- now reset the coeffarray for the next allocation
         coeffarray(icoeff) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_tornames(itype,nsame,namearray)
c     ******************************************************************
c     * assigns all of the tornames and sets ntsame for this type of   *
c     * torsion angle                                                  *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMET,MAXTNAME)
c     --- local variables
      integer isame,iatom

      call twh_ntsame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do iatom = 1,MAXTNAME
            call twh_tornames(GLB_SET,itype,isame,iatom
     &           ,namearray(isame,iatom))
c           --- now reset the namearray for the next allocation
            namearray(isame,iatom) = '-'
         enddo
      enddo

      return
      end

      subroutine twh_set_vibcoeffs(itype,vibarray)
c     ******************************************************************
c     * assigns all of the vibcoeffs for a single vibration type       *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision vibarray
      dimension vibarray(MINVIBCOEFF:MAXVIBCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINVIBCOEFF,MAXVIBCOEFF
         call twh_vibcoeff(GLB_SET,itype,icoeff,vibarray(icoeff))
      enddo

      return
      end

      subroutine twh_set_vibnames(itype,nsame,namearray)
c     ******************************************************************
c     * assigns all of the vibnames and sets nvsame for this type of   *
c     * vibration                                                      *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMEV,MAXVNAME)
c     --- local variables
      integer isame,iatom

      call twh_nvsame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do iatom = 1,MAXVNAME
            call twh_vibnames(GLB_SET,itype,isame,iatom
     &           ,namearray(isame,iatom))
c           --- now reset the namearray for the next allocation
            namearray(isame,iatom) = '-'
         enddo
      enddo

      return
      end

c     ******************************************************************
c     * functions that primarily make it easier to access the global   *
c     * data structures in a manner that is similar to accessing an    *
c     * array                                                          *
c     * also includes functions that resemble the fetching of global   *
c     * data structures                                                *
c     ******************************************************************

      function twh_get_aacoeff(itype,icoeff)
c     ******************************************************************
c     * gets the aacoeff data                                          *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_AACOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,icoeff

      call twh_aacoeff(GLB_GET,itype,icoeff,twh_get_aacoeff)
  
      return
      end

      function twh_get_aanames(itype,isame,index)
c     ******************************************************************
c     * gets the aanames data                                          *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_AANAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,isame,index

      call twh_aanames(GLB_GET,itype,isame,index,twh_get_aanames)
  
      return
      end

      function twh_get_aaff(itype)
c     ******************************************************************
c     * gets the aaff data                                             *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_AAFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_aaff(GLB_GET,itype,twh_get_aaff)
  
      return
      end

      function twh_get_aastyle(itype)
c     ******************************************************************
c     * gets the aastyle data                                          *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_AASTYLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_aastyle(GLB_GET,itype,twh_get_aastyle)
  
      return
      end

      function twh_get_anglestyle(itype)
c     ******************************************************************
c     * gets the anglestyle data                                       *
c     *                                                                *
c     * originally written 04-04-2006 by M.G. Martin                   *
c     * last modified 04-04-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ANGLESTYLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_anglestyle(GLB_GET,itype,twh_get_anglestyle)
  
      return
      end

      function twh_get_bencoeff(type,icoeff)
c     ******************************************************************
c     * gets the bencoeff data                                         *
c     *                                                                *
c     * originally written 04-04-2006 by M.G. Martin                   *
c     * last modified 04-04-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_BENCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer type,icoeff

      call twh_bencoeff(GLB_GET,type,icoeff,twh_get_bencoeff)
  
      return
      end

      function twh_get_bendff(type)
c     ******************************************************************
c     * gets the bendff data                                           *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_BENDFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer type

      call twh_bendff(GLB_GET,type,twh_get_bendff)
  
      return
      end

      function twh_get_bendnames(type,isame,index)
c     ******************************************************************
c     * gets the bendnames data                                        *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_BENDNAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer type,isame,index

      call twh_bendnames(GLB_GET,type,isame,index,twh_get_bendnames)
  
      return
      end

      function twh_get_benprefact(imolty,iunit,ibend)
c     ******************************************************************
c     * gets the benprefact data                                       *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_BENPREFACT
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibend

      call twh_benprefact(GLB_GET,imolty,iunit,ibend,twh_get_benprefact)
  
      return
      end

      function twh_get_bondpatt(itype)
c     ******************************************************************
c     * gets the bondpatt data                                         *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_BONDPATT
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_bondpatt(GLB_GET,itype,twh_get_bondpatt)
  
      return
      end

      function twh_get_caname(itype)
c     ******************************************************************
c     * gets the caname data                                           *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CANAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_caname(GLB_GET,itype,twh_get_caname)
  
      return
      end

      function twh_get_cbname(itype)
c     ******************************************************************
c     * gets the cbname data                                           *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CBNAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_cbname(GLB_GET,itype,twh_get_cbname)
  
      return
      end

      function twh_get_classical_mixrule()
c     ******************************************************************
c     * gets the classical_mixrule data                                *
c     *                                                                *
c     * originally written 04-15-2006 by M.G. Martin                   *
c     * last modified 04-15-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CLASSICAL_MIXRULE
#include "functions.h"

      call twh_classical_mixrule(GLB_GET,twh_get_classical_mixrule)
  
      return
      end

      function twh_get_classical_potential()
c     ******************************************************************
c     * gets the classical_potential data                              *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#include "functions.h"

      call twh_classical_potential(GLB_GET,twh_get_classical_potential)
  
      return
      end

      function twh_get_ctname(itype)
c     ******************************************************************
c     * gets the ctname data                                           *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CTNAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_ctname(GLB_GET,itype,twh_get_ctname)
  
      return
      end

      function twh_get_elementname(itype)
c     ******************************************************************
c     * gets the elementname data                                      *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ELEMENTNAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_elementname(GLB_GET,itype,twh_get_elementname)
  
      return
      end

      function twh_get_growcount(imolty,iunit,ibond)
c     ******************************************************************
c     * gets the growcount data                                        *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_GROWCOUNT
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibond

      call twh_growcount(GLB_GET,imolty,iunit,ibond,twh_get_growcount)
  
      return
      end

      function twh_get_growvalidcount(imolty)
c     ******************************************************************
c     * gets the growvalidcount data                                   *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_GROWVALIDCOUNT
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty

      call twh_growvalidcount(GLB_GET,imolty,twh_get_growvalidcount)
  
      return
      end

      function twh_get_growvalidlist(imolty,icount)
c     ******************************************************************
c     * gets the growvalidlist data                                    *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_GROWVALIDLIST
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,icount

      call twh_growvalidlist(GLB_GET,imolty,icount
     &     ,twh_get_growvalidlist)
  
      return
      end

      function twh_get_hbondcoeff(itype,index)
c     ******************************************************************
c     * gets the hbondcoeff data                                       *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_HBONDCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,index

      call twh_hbondcoeff(GLB_GET,itype,index,twh_get_hbondcoeff)
  
      return
      end

      function twh_get_hbondnames(itype,index)
c     ******************************************************************
c     * gets the hbondnames data                                       *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_HBONDNAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,index

      call twh_hbondnames(GLB_GET,itype,index,twh_get_hbondnames)
  
      return
      end

      function twh_get_ibtben1(imolty,iunit,ibend)
c     ******************************************************************
c     * gets the ibtben1 data                                          *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IBTBEN1
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibend

      call twh_ibtben1(GLB_GET,imolty,iunit,ibend,twh_get_ibtben1)
  
      return
      end

      function twh_get_ibtben2(imolty,iunit,ibend)
c     ******************************************************************
c     * gets the ibtben2 data                                          *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IBTBEN2
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibend

      call twh_ibtben2(GLB_GET,imolty,iunit,ibend,twh_get_ibtben2)
  
      return
      end

      function twh_get_ijaa0(imolty,iunit,iaa)
c     ******************************************************************
c     * gets the ijaa0 data                                            *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJAA0
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iaa

      call twh_ijaa0(GLB_GET,imolty,iunit,iaa,twh_get_ijaa0)
  
      return
      end

      function twh_get_ijaa1(imolty,iunit,iaa)
c     ******************************************************************
c     * gets the ijaa1 data                                            *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJAA1
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iaa

      call twh_ijaa1(GLB_GET,imolty,iunit,iaa,twh_get_ijaa1)
  
      return
      end

      function twh_get_ijaa2(imolty,iunit,iaa)
c     ******************************************************************
c     * gets the ijaa2 data                                            *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJAA2
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iaa

      call twh_ijaa2(GLB_GET,imolty,iunit,iaa,twh_get_ijaa2)
  
      return
      end

      function twh_get_ijben2(imolty,iunit,ibend)
c     ******************************************************************
c     * gets the ijben2 data                                           *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJBEN2
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibend

      call twh_ijben2(GLB_GET,imolty,iunit,ibend,twh_get_ijben2)
  
      return
      end

      function twh_get_ijben3(imolty,iunit,ibend)
c     ******************************************************************
c     * gets the ijben3 data                                           *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJBEN3
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibend

      call twh_ijben3(GLB_GET,imolty,iunit,ibend,twh_get_ijben3)
  
      return
      end

      function twh_get_ijbnbn1(imolty,iunit,ibend)
c     ******************************************************************
c     * gets the ijbnbn1 data                                          *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJBNBN1
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibend

      call twh_ijbnbn1(GLB_GET,imolty,iunit,ibend,twh_get_ijbnbn1)
  
      return
      end

      function twh_get_ijbnbn2(imolty,iunit,ibend)
c     ******************************************************************
c     * gets the ijbnbn2 data                                          *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJBNBN2
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibend

      call twh_ijbnbn2(GLB_GET,imolty,iunit,ibend,twh_get_ijbnbn2)
  
      return
      end

      function twh_get_ijimprop2(imolty,iunit,iimp)
c     ******************************************************************
c     * gets the ijimprop2 data                                        *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJIMPROP2
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iimp

      call twh_ijimprop2(GLB_GET,imolty,iunit,iimp,twh_get_ijimprop2)
  
      return
      end

      function twh_get_ijimprop3(imolty,iunit,iimp)
c     ******************************************************************
c     * gets the ijimprop3 data                                        *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJIMPROP3
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iimp

      call twh_ijimprop3(GLB_GET,imolty,iunit,iimp,twh_get_ijimprop3)
  
      return
      end

      function twh_get_ijimprop4(imolty,iunit,iimp)
c     ******************************************************************
c     * gets the ijimprop4 data                                        *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJIMPROP4
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iimp

      call twh_ijimprop4(GLB_GET,imolty,iunit,iimp,twh_get_ijimprop4)
  
      return
      end

      function twh_get_ijof5(imolty,iunit,iof)
c     ******************************************************************
c     * gets the ijof5 data                                            *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJOF5
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iof

      call twh_ijof5(GLB_GET,imolty,iunit,iof,twh_get_ijof5)
  
      return
      end

      function twh_get_ijtor2(imolty,iunit,itor)
c     ******************************************************************
c     * gets the ijtor2 data                                           *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJTOR2
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,itor

      call twh_ijtor2(GLB_GET,imolty,iunit,itor,twh_get_ijtor2)
  
      return
      end

      function twh_get_ijtor3(imolty,iunit,itor)
c     ******************************************************************
c     * gets the ijtor3 data                                           *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJTOR3
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,itor

      call twh_ijtor3(GLB_GET,imolty,iunit,itor,twh_get_ijtor3)
  
      return
      end

      function twh_get_ijtor4(imolty,iunit,itor)
c     ******************************************************************
c     * gets the ijtor4 data                                           *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJTOR4
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,itor

      call twh_ijtor4(GLB_GET,imolty,iunit,itor,twh_get_ijtor4)
  
      return
      end

      function twh_get_ijvib(imolty,iunit,ibond)
c     ******************************************************************
c     * gets the ijvib data                                            *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJVIB
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibond

      call twh_ijvib(GLB_GET,imolty,iunit,ibond,twh_get_ijvib)
  
      return
      end

      function twh_get_impcoeff(itype,icoeff)
c     ******************************************************************
c     * gets the impcoeff data                                         *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,icoeff

      call twh_impcoeff(GLB_GET,itype,icoeff,twh_get_impcoeff)
  
      return
      end

      function twh_get_impff(itype)
c     ******************************************************************
c     * gets the impff data                                            *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_impff(GLB_GET,itype,twh_get_impff)
  
      return
      end

      function twh_get_impform(itype)
c     ******************************************************************
c     * gets the impform data                                          *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPFORM
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_impform(GLB_GET,itype,twh_get_impform)
  
      return
      end

      function twh_get_impnames(itype,isame,index)
c     ******************************************************************
c     * gets the impnames data                                         *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPNAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,isame,index

      call twh_impnames(GLB_GET,itype,isame,index,twh_get_impnames)
  
      return
      end

      function twh_get_impstyle(itype)
c     ******************************************************************
c     * gets the impstyle data                                         *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPSTYLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_impstyle(GLB_GET,itype,twh_get_impstyle)
  
      return
      end

      function twh_get_inaa(imolty,iunit)
c     ******************************************************************
c     * gets the inaa data                                             *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INAA
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit

      call twh_inaa(GLB_GET,imolty,iunit,twh_get_inaa)
  
      return
      end

      function twh_get_inben(imolty,iunit)
c     ******************************************************************
c     * gets the inben data                                            *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INBEN
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit

      call twh_inben(GLB_GET,imolty,iunit,twh_get_inben)
  
      return
      end

      function twh_get_inbnbn(imolty,iunit)
c     ******************************************************************
c     * gets the inbnbn data                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INBNBN
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit

      call twh_inbnbn(GLB_GET,imolty,iunit,twh_get_inbnbn)
  
      return
      end

      function twh_get_inimprop(imolty,iunit)
c     ******************************************************************
c     * gets the inimprop data                                         *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INIMPROP
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit

      call twh_inimprop(GLB_GET,imolty,iunit,twh_get_inimprop)
  
      return
      end

      function twh_get_inof(imolty,iunit)
c     ******************************************************************
c     * gets the inof data                                             *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INOF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit

      call twh_inof(GLB_GET,imolty,iunit,twh_get_inof)
  
      return
      end

      function twh_get_intor(imolty,iunit)
c     ******************************************************************
c     * gets the intor data                                            *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INTOR
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit

      call twh_intor(GLB_GET,imolty,iunit,twh_get_intor)
  
      return
      end

      function twh_get_invib(imolty,iunit)
c     ******************************************************************
c     * gets the invib data                                            *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INVIB
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit

      call twh_invib(GLB_GET,imolty,iunit,twh_get_invib)
  
      return
      end

      function twh_get_itaa(imolty,iunit,iaa)
c     ******************************************************************
c     * gets the itaa data                                             *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ITAA
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iaa

      call twh_itaa(GLB_GET,imolty,iunit,iaa,twh_get_itaa)
  
      return
      end

      function twh_get_itben(imolty,iunit,ibend)
c     ******************************************************************
c     * gets the itben data                                            *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ITBEN
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibend

      call twh_itben(GLB_GET,imolty,iunit,ibend,twh_get_itben)
  
      return
      end

      function twh_get_itbnbn(imolty,iunit,ibend)
c     ******************************************************************
c     * gets the itbnbn data                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ITBNBN
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibend

      call twh_itbnbn(GLB_GET,imolty,iunit,ibend,twh_get_itbnbn)
  
      return
      end

      function twh_get_itimprop(imolty,iunit,iimp)
c     ******************************************************************
c     * gets the itimprop data                                         *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ITIMPROP
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iimp

      call twh_itimprop(GLB_GET,imolty,iunit,iimp,twh_get_itimprop)
  
      return
      end

      function twh_get_itof(imolty,iunit,iof)
c     ******************************************************************
c     * gets the itof data                                             *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ITOF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,iof

      call twh_itof(GLB_GET,imolty,iunit,iof,twh_get_itof)
  
      return
      end

      function twh_get_ittor(imolty,iunit,itor)
c     ******************************************************************
c     * gets the ittor data                                            *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ITTOR
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,itor

      call twh_ittor(GLB_GET,imolty,iunit,itor,twh_get_ittor)
  
      return
      end

      function twh_get_itvib(imolty,iunit,ibond)
c     ******************************************************************
c     * gets the itvib data                                            *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ITVIB
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iunit,ibond

      call twh_itvib(GLB_GET,imolty,iunit,ibond,twh_get_itvib)
  
      return
      end

      function twh_get_laahere(iaatype)
c     ******************************************************************
c     * returns true if this angle-angle type appears in a molecule    *
c     * that is being used in this simulation                          *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INAA
#define FUNCTION_GET_ITAA
#define FUNCTION_GET_LAAHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer iaatype
c     --- local variables
      integer imolty,iunit,iaa,nmolty,nunit

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            do iaa = 1,twh_get_inaa(imolty,iunit)
               if ( twh_get_itaa(imolty,iunit,iaa) .eq. iaatype ) then
                  twh_get_laahere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_laahere = .false.

      return
      end

      function twh_get_lbendhere(bendtype)
c     ******************************************************************
c     * returns true if this angle type appears in a molecule that     *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LBENDHERE
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_ITBEN
#include "functions.h"
c     --- variables passed to/from the function
      integer bendtype
c     --- local variables
      integer imolty,iunit,ibend,nmolty,nunit

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            do ibend = 1,twh_get_inben(imolty,iunit)
               if ( twh_get_itben(imolty,iunit,ibend)
     &              .eq. bendtype ) then
                  twh_get_lbendhere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_lbendhere = .false.

      return
      end

      function twh_get_lbonang(ibtype)
c     ******************************************************************
c     * gets the lbonang data                                          *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LBONANG
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ibtype

      call twh_lbonang(GLB_GET,ibtype,twh_get_lbonang)
  
      return
      end

      function twh_get_lbonbon(ibtype)
c     ******************************************************************
c     * gets the lbonbon data                                          *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LBONBON
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ibtype

      call twh_lbonbon(GLB_GET,ibtype,twh_get_lbonbon)
  
      return
      end

      function twh_get_lexist(iunit)
c     ******************************************************************
c     * gets the lexist data                                           *
c     *                                                                *
c     * originally written 12-14-2007 by M.G. Martin                   *
c     * last modified 12-14-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LEXIST
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer iunit

      call twh_lexist(GLB_GET,iunit,twh_get_lexist)
  
      return
      end

      function twh_get_lhere(nbtype)
c     ******************************************************************
c     * returns true if this nonbonded type appears in a molecule that *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer nbtype
c     --- local variables
      integer imolty,iunit,nmolty,nunit,itype

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_ntype(GLB_GET,imolty,iunit,itype)
            if ( itype .eq. nbtype ) then
               twh_get_lhere = .true.
               return
            endif
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_lhere = .false.

      return
      end

      function twh_get_limphere(imptype)
c     ******************************************************************
c     * returns true if this improper type appears in a molecule that  *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LIMPHERE
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_ITIMPROP
#include "functions.h"
c     --- variables passed to/from the function
      integer imptype
c     --- local variables
      integer imolty,iunit,iimp,nmolty,nunit

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            do iimp = 1,twh_get_inimprop(imolty,iunit)
               if ( twh_get_itimprop(imolty,iunit,iimp) 
     &              .eq. imptype ) then
                  twh_get_limphere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_limphere = .false.

      return
      end

      function twh_get_lofhere(oftype)
c     ******************************************************************
c     * returns true if this one-five type appears in a molecule       *
c     * that is being used in this simulation                          *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INOF
#define FUNCTION_GET_ITOF
#define FUNCTION_GET_LOFHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer oftype
c     --- local variables
      integer imolty,iunit,iof,nmolty,nunit

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            do iof = 1,twh_get_inof(imolty,iunit)
               if ( twh_get_itof(imolty,iunit,iof) .eq. oftype ) then
                  twh_get_lofhere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_lofhere = .false.

      return
      end

      function twh_get_loftor(ittype)
c     ******************************************************************
c     * gets the loftor data                                           *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LOFTOR
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ittype

      call twh_loftor(GLB_GET,ittype,twh_get_loftor)
  
      return
      end

      function twh_get_ltorhere(tortype)
c     ******************************************************************
c     * returns true if this torsion type appears in a molecule that   *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LTORHERE
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_ITTOR
#include "functions.h"
c     --- variables passed to/from the function
      integer tortype
c     --- local variables
      integer imolty,iunit,itor,nmolty,nunit

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            do itor = 1,twh_get_intor(imolty,iunit)
               if ( abs(twh_get_ittor(imolty,iunit,itor)) 
     &              .eq. abs(tortype) ) then
                  twh_get_ltorhere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_ltorhere = .false.

      return
      end

      function twh_get_lvibhere(vibtype)
c     ******************************************************************
c     * returns true if this vibration type appears in a molecule that *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LVIBHERE
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_ITVIB
#include "functions.h"
c     --- variables passed to/from the function
      integer vibtype
c     --- local variables
      integer imolty,iunit,ivib,nmolty,nunit

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            do ivib = 1,twh_get_invib(imolty,iunit)
               if ( twh_get_itvib(imolty,iunit,ivib) .eq. vibtype ) then
                  twh_get_lvibhere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_lvibhere = .false.

      return
      end

      function twh_get_mass(itype)
c     ******************************************************************
c     * gets the mass data                                             *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_MASS
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_mass(GLB_GET,itype,twh_get_mass)
  
      return
      end

      function twh_get_naasame(itype)
c     ******************************************************************
c     * gets the naasame data                                          *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NAASAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_naasame(GLB_GET,itype,twh_get_naasame)
  
      return
      end

      function twh_get_nbcoeff(itype,icoeff)
c     ******************************************************************
c     * gets the nbcoeff data                                          *
c     *                                                                *
c     * originally written 04-11-2006 by M.G. Martin                   *
c     * last modified 04-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,icoeff

      call twh_nbcoeff(GLB_GET,itype,icoeff,twh_get_nbcoeff)
  
      return
      end

      function twh_get_nbff(itype)
c     ******************************************************************
c     * gets the nbff data                                             *
c     *                                                                *
c     * originally written 04-11-2006 by M.G. Martin                   *
c     * last modified 04-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_nbff(GLB_GET,itype,twh_get_nbff)
  
      return
      end

      function twh_get_nbname(itype)
c     ******************************************************************
c     * gets the nbname data                                           *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBNAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_nbname(GLB_GET,itype,twh_get_nbname)
  
      return
      end

      function twh_get_nboxi(ichain)
c     ******************************************************************
c     * gets the nboxi data                                            *
c     *                                                                *
c     * originally written 11-01-2007 by M.G. Martin                   *
c     * last modified 11-01-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBOXI
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain

      call twh_nboxi(GLB_GET,ichain,twh_get_nboxi)
  
      return
      end

      function twh_get_nbsame(itype)
c     ******************************************************************
c     * gets the nbsame data                                           *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBSAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_nbsame(GLB_GET,itype,twh_get_nbsame)
  
      return
      end

      function twh_get_ncmt(ibox,itype)
c     ******************************************************************
c     * gets the ncmt data                                             *
c     *                                                                *
c     * originally written 11-01-2007 by M.G. Martin                   *
c     * last modified 11-01-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NCMT
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ibox,itype

      call twh_ncmt(GLB_GET,ibox,itype,twh_get_ncmt)
  
      return
      end

      function twh_get_nimpsame(itype)
c     ******************************************************************
c     * gets the nimpsame data                                         *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NIMPSAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_nimpsame(GLB_GET,itype,twh_get_nimpsame)
  
      return
      end

      function twh_get_ntsame(itype)
c     ******************************************************************
c     * gets the ntsame data                                           *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NTSAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_ntsame(GLB_GET,itype,twh_get_ntsame)
  
      return
      end

      function twh_get_ntorloop(itype)
c     ******************************************************************
c     * gets the ntorloop data                                         *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NTORLOOP
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_ntorloop(GLB_GET,itype,twh_get_ntorloop)
  
      return
      end

      function twh_get_ntype(itype,iunit)
c     ******************************************************************
c     * gets the ntype data                                            *
c     *                                                                *
c     * originally written 11-01-2007 by M.G. Martin                   *
c     * last modified 11-01-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,iunit

      call twh_ntype(GLB_GET,itype,iunit,twh_get_ntype)
  
      return
      end

      function twh_get_nvsame(itype)
c     ******************************************************************
c     * gets the nvsame data                                           *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NVSAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_nvsame(GLB_GET,itype,twh_get_nvsame)
  
      return
      end

      function twh_get_ofcoeff(itype,icoeff)
c     ******************************************************************
c     * gets the ofcoeff data                                          *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_OFCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,icoeff

      call twh_ofcoeff(GLB_GET,itype,icoeff,twh_get_ofcoeff)
  
      return
      end

      function twh_get_offf(itype)
c     ******************************************************************
c     * gets the offf data                                             *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_OFFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_offf(GLB_GET,itype,twh_get_offf)
  
      return
      end

      function twh_get_ofnames(itype,index)
c     ******************************************************************
c     * gets the ofnames data                                          *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_OFNAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,index

      call twh_ofnames(GLB_GET,itype,index,twh_get_ofnames)
  
      return
      end

      function twh_get_ofstyle(itype)
c     ******************************************************************
c     * gets the ofstyle data                                          *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_OFSTYLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_ofstyle(GLB_GET,itype,twh_get_ofstyle)
  
      return
      end

      function twh_get_pm1boxcomswitch()
c     ******************************************************************
c     * gets the pm1boxcomswitch data                                  *
c     *                                                                *
c     * originally written 07-28-2006 by M.G. Martin                   *
c     * last modified 07-28-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_PM1BOXCOMSWITCH
#include "functions.h"

      call twh_pm1boxcomswitch(GLB_GET,twh_get_pm1boxcomswitch)
  
      return
      end

      function twh_get_pmvlpr(ipairbox)
c     ******************************************************************
c     * gets the pmvlpr data                                           *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_PMVLPR
#include "functions.h"
c     --- variables passed to/from the function
      integer ipairbox

      call twh_pmvlpr(GLB_GET,ipairbox,twh_get_pmvlpr)
  
      return
      end

      function twh_get_pmvol()
c     ******************************************************************
c     * gets the pmvol data                                            *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_PMVOL
#include "functions.h"

      call twh_pmvol(GLB_GET,twh_get_pmvol)
  
      return
      end

      function twh_get_qbasevalue(itype)
c     ******************************************************************
c     * gets the qbasevalue data                                       *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_QBASEVALUE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_qbasevalue(GLB_GET,itype,twh_get_qbasevalue)
  
      return
      end

      function twh_get_qbiff(itype)
c     ******************************************************************
c     * gets the qbiff data                                            *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_QBIFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_qbiff(GLB_GET,itype,twh_get_qbiff)
  
      return
      end

      function twh_get_qbinames(itype,index)
c     ******************************************************************
c     * gets the qbinames data                                         *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_QBINAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,index

      call twh_qbinames(GLB_GET,itype,index,twh_get_qbinames)
  
      return
      end

      function twh_get_qbivalue(itype)
c     ******************************************************************
c     * gets the qbivalue data                                         *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_QBIVALUE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_qbivalue(GLB_GET,itype,twh_get_qbivalue)
  
      return
      end

      function twh_get_qqatom(itype,iunit)
c     ******************************************************************
c     * gets the qqatom data                                           *
c     *                                                                *
c     * originally written 11-02-2007 by M.G. Martin                   *
c     * last modified 11-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_QQATOM
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,iunit

      call twh_qqatom(GLB_GET,itype,iunit,twh_get_qqatom)
  
      return
      end

      function twh_get_rmvol(ipairbox)
c     ******************************************************************
c     * gets the rmvol data                                            *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_RMVOL
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ipairbox

      call twh_rmvol(GLB_GET,ipairbox,twh_get_rmvol)
  
      return
      end

      function twh_get_tavol()
c     ******************************************************************
c     * gets the tavol data                                            *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_TAVOL
#include "functions.h"

      call twh_tavol(GLB_GET,twh_get_tavol)
  
      return
      end

      function twh_get_torff(itype)
c     ******************************************************************
c     * gets the torff data                                            *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_TORFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_torff(GLB_GET,itype,twh_get_torff)
  
      return
      end

      function twh_get_torcoeff(itype,icoeff)
c     ******************************************************************
c     * gets the torcoeff data                                         *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_TORCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,icoeff

      call twh_torcoeff(GLB_GET,itype,icoeff,twh_get_torcoeff)
  
      return
      end

      function twh_get_tornames(itype,isame,index)
c     ******************************************************************
c     * gets the tornames data                                         *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_TORNAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,isame,index

      call twh_tornames(GLB_GET,itype,isame,index,twh_get_tornames)
  
      return
      end

      function twh_get_torstring(itype,index)
c     ******************************************************************
c     * gets the torstring data                                        *
c     *                                                                *
c     * originally written 04-06-2006 by M.G. Martin                   *
c     * last modified 04-06-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_TORSTRING
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,index

      call twh_torstring(GLB_GET,itype,index,twh_get_torstring)
  
      return
      end

      function twh_get_torstyle(itype)
c     ******************************************************************
c     * gets the torstyle data                                         *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_TORSTYLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype

      call twh_torstyle(GLB_GET,itype,twh_get_torstyle)
  
      return
      end

      function twh_get_vibcoeff(type,icoeff)
c     ******************************************************************
c     * gets the vibcoeff data                                         *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_VIBCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer type,icoeff

      call twh_vibcoeff(GLB_GET,type,icoeff,twh_get_vibcoeff)
  
      return
      end

      function twh_get_vibff(type)
c     ******************************************************************
c     * gets the vibff data                                            *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_VIBFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer type

      call twh_vibff(GLB_GET,type,twh_get_vibff)
  
      return
      end

      function twh_get_vibnames(itype,isame,iatom)
c     ******************************************************************
c     * gets the vibnames data                                         *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_VIBNAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,isame,iatom

      call twh_vibnames(GLB_GET,itype,isame,iatom,twh_get_vibnames)
  
      return
      end

      function twh_get_vibrang(index)
c     ******************************************************************
c     * gets the vibrang data                                          *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_VIBRANG
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer index

      call twh_vibrang(GLB_GET,index,twh_get_vibrang)
  
      return
      end


      function twh_limited_double(dvalue,limit)
c     ******************************************************************
c     * check that the magnitude of a parameters has an absolute value *
c     * that is sufficiently large to reduce the noise between machines*
c     *                                                                *
c     * originally written 12-17-2007 by M.G. Martin                   *
c     * last modified 12-22-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_LIMITED_DOUBLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision dvalue,limit

      if ( dabs(dvalue) .lt. limit ) then
         twh_limited_double = 0.0d0
      else
         twh_limited_double = dvalue
      endif
  
      return
      end

      function twh_safe_double(dvalue)
c     ******************************************************************
c     * checks to see that the double precision dvalue is larger than  *
c     * SMALLEST.  If it is larger then it returns the value of dvalue *
c     * If it is smaller than SMALLEST it returns 0.0d0                *
c     * This is used when outputing information to reduce machine      *
c     * precision noise                                                *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_SAFE_DOUBLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision dvalue

      if ( dabs(dvalue) .lt. SMALLEST ) then
         twh_safe_double = 0.0d0
      else
         twh_safe_double = dvalue
      endif
  
      return
      end

