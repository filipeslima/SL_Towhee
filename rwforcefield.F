#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C)2002-2014 Marcus G. Martin                        *
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
      subroutine twh_rwforcefield(lfailure,rwflag,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes,nimprotypes,naatypes,nbitypes
     &     ,noftypes,nhbondtypes,iffnum)
c     ******************************************************************
c     * reads or writes all of the force field information to/from     *
c     * towhee_ff in a highly formatted fashion                        *
c     * rwflag must be either READ_FLAG or WRITE_FLAG                  *
c     *                                                                *
c     * rewritten from readforcefield and writeforcefield 04-01-2002   *
c     * by M.G. Martin                                                 *
c     * last modified 07-18-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CLASSICAL_MIXRULE
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#define FUNCTION_GET_NCOEFF
#define FUNCTION_READ_INTEGER
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_READ_LABELED_LOGICAL
#include "functions.h"
c     --- local definitions
#define RWFORCEFIELD_DEBUG 0 
#define FFUNIT 17
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
      integer naatypes,iffnum,nbitypes,noftypes,nhbondtypes
c     --- local variables
      logical lexplicit,lgetcoeffnum
      logical lrwhbond
      integer type, currentffversion,iversion,ncoeff
      integer nb_adj,bond_adj,angle_adj,dihed_adj,impr_adj,aa_adj,bi_adj
      integer of_adj,hb_adj
      integer hbcoeffstart,hbcoeffstop
      character*15 idname
      character*30 localptype,oldmtype,cmlocal
      character*50 goodstring
      character*(MAXDIRLENGTH) fffilename

      if (rwflag .ne. WRITE_FLAG .and. rwflag .ne. READ_FLAG) then
         write(6,*) "Illegal rwflag in rwforcefield"
         stop
      endif
#if RWFORCEFIELD_DEBUG
      if (rwflag .eq. WRITE_FLAG)
         write(6,*) 'starting RWFORCEFIELD -- writing'
      else
         write(6,*) 'starting RWFORCEFIELD -- reading'
#endif

c     --- set subroutine idname
      idname = 'rwforcefield'
c     --- set the current version of towhee_ff
      currentffversion = 15
c     --- open the forcefield file
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- writing a forcefield file
c        --- prepend with the directory in the usual manner
         call twh_ff_filename(GLB_GET,iffnum,fffilename)
         call twh_openfile(lfailure,FFUNIT,fffilename,'DEFAULT'
     &        ,'UNKNOWN','DEFAULT')
         if ( lfailure ) return
      else
c        --- reading a forcefield file
c        --- this is an exception to the general rule that prepends
c        --- the directory to every file.
         call twh_ff_filename(GLB_GET,iffnum,fffilename)
         open (unit=FFUNIT,file=fffilename,form='formatted',
     &      status='OLD',ERR=100)
      endif

c     --- set the adjustments
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- we are writing to a single file.  all adjustments are 0
         nb_adj = 0
         bond_adj = 0
         bi_adj = 0
         angle_adj = 0
         dihed_adj = 0
         impr_adj = 0
         aa_adj = 0
         of_adj = 0
         hb_adj = 0
      else
c        --- could be reading multiple files, adjustments are passed in
         nb_adj = ntypes
         bond_adj = nbondtypes
         bi_adj = nbitypes
         angle_adj = nangletypes
         dihed_adj = ndihedtypes
         impr_adj = nimprotypes
         aa_adj = naatypes
         of_adj = noftypes
         hb_adj = nhbondtypes
      endif

c     --- read the towhee_ff version
      goodstring = 'towhee_ff Version'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a17)') goodstring
         write(FFUNIT,*) currentffversion
c        --- set iversion to the currentffversion
         iversion = currentffversion
      else
         iversion =  twh_read_labeled_integer(lfailure,FFUNIT
     &        ,goodstring,.false.,idname)
         if ( lfailure ) return

c        --- cannot work with version earlier than 4
         if ( iversion .lt. 4 ) then
            write(6,*) 'RWFORCEFIELD: No longer compatable with '
     &           ,'towhee_ff versions less than 4'
            stop
         endif
         if ( iversion .gt. currentffversion ) then
            write(6,*) 'invalid towhee_ff version ',iversion
            stop
         endif
      endif
      
c     --- number of nonbonded types
      goodstring = 'Number of Nonbonded Types'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a25)') goodstring
         write(FFUNIT,*) ntypes
      else
         ntypes =  twh_read_labeled_integer(lfailure,FFUNIT,goodstring,
     &         .false.,idname)
         if ( lfailure ) return
         if ( ntypes + nb_adj .gt. NNTYPE ) then
            write(6,*) 'ntypes,nb_adj, NNTYPE ',ntypes,nb_adj,NNTYPE
            write(6,*) 'RWFORCEFIELD: Number of Nonbonded Types '
     &           ,'exceeds NNTYPE'
            stop
         endif
      endif

c     --- the nonbonded potential
      if ( rwflag .eq. WRITE_FLAG ) then
         goodstring = 'Potential Type'
         write(FFUNIT,'(a14)') goodstring
         write(FFUNIT,'(a30)') twh_get_classical_potential()
c        --- set the localptype to the classical_potential
         localptype = twh_get_classical_potential()
      else
         call twh_get_potential_type(lfailure,iversion, localptype)
         if ( lfailure ) return
      endif

      goodstring = 'Classical Mixrule'
      call twh_classical_mixrule(GLB_GET,cmlocal)
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a17)') goodstring
         write(FFUNIT,'(a30)') cmlocal
      else
         if ( iversion .ge. 7 ) then
c           --- read in the mixrule and error check
            call twh_read_labeled_string(lfailure,FFUNIT,goodstring,
     &           .false.,idname,oldmtype)
            if ( lfailure ) return
c           --- error check oldmtype against classical_mixrule
            if ( oldmtype .ne. cmlocal 
     &           .and.
     &           .not. ( oldmtype .eq. 'LB or Geometric' .and.
     &           ( cmlocal .eq. 'Lorentz-Berthelot' .or.
     &           cmlocal .eq. 'Geometric' .or.
     &           cmlocal .eq. 'LB plus manual' ) )
     &           .and.
     &           .not. ( oldmtype .eq. 'Lorentz-Berthelot' .and.
     &           cmlocal .eq. 'LB plus manual' )
     &           .and. 
     &           ( oldmtype .ne. 'wild' )
     &           ) then
               write(6,*) 'RWFORCEFIELD: mismatch of classical_mixrule'
               write(6,*) 'towhee_input: ',cmlocal
               write(6,*) 'ff_filename: ',oldmtype
               stop
            endif
         endif
      endif
      
c     --- for many of the potentials the number of nbcoeffs is constant
c     --- in that case, this number is set here.  Otherwise it is set
c     --- later
      lgetcoeffnum = .false.
c     --- initialze the hbond logical to false as it is almost never 
c     --- used.  If it is needed it will be reset to true in get_ncoeff
      lrwhbond = .false.

      ncoeff = twh_get_ncoeff(localptype, iversion, lrwhbond
     &     , lgetcoeffnum, hbcoeffstart, hbcoeffstop)

      if ( iversion .lt. 7 ) then
c        --- explicit cross terms logical
         goodstring = 'Cross Term Logical'
         if ( rwflag .eq. WRITE_FLAG ) then
            if ( twh_get_classical_mixrule() .eq. 'Explicit' ) then
               lexplicit = .true.
            else
               lexplicit = .false.
            endif
            write(FFUNIT,'(a18)') goodstring
            write(FFUNIT,*) lexplicit
         else
            lexplicit =  twh_read_labeled_logical(lfailure,FFUNIT
     &           ,goodstring, .false.,idname)
            if ( lfailure ) return
         endif
      else
c        --- set cross term logical based upon classical_mixrule
         if ( twh_get_classical_mixrule() .eq. 'Explicit' ) then
            lexplicit = .true.
         else
            lexplicit = .false.
         endif
      endif

      do type = nb_adj+1,nb_adj+ntypes
         call twh_rw_atom_data(lfailure,rwflag, type, nb_adj, ntypes
     &        ,iffnum,iversion,lexplicit,lgetcoeffnum,localptype,ncoeff)
         if ( lfailure ) return
      enddo

c     --- bond force field information
      goodstring = 'Number of Bonded Terms'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a22)') goodstring
         write(FFUNIT,*) nbondtypes
      else
         nbondtypes =  twh_read_labeled_integer(lfailure, FFUNIT,
     &         goodstring, .false.,idname)
         if ( lfailure ) return
      endif
      if ( nbondtypes + bond_adj .gt. TVIBMAX ) then
         write(6,*) 'RWFORCEFIELD: Number of Bonded Terms exceeds '
     &        ,' TVIBMAX'
         write(6,*) 'nbondtypes,bond_adj TVIBMAX',nbondtypes
     &        ,bond_adj,TVIBMAX
         stop
      endif

      do type = bond_adj+1,bond_adj+nbondtypes
         call twh_rw_bond_data(lfailure,iversion,rwflag,type,bond_adj)
         if ( lfailure ) return
      enddo

c     --- Angle force field information
      goodstring = 'Number of Angle Terms'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a21)') goodstring
         write(FFUNIT,*) nangletypes
      else
         nangletypes =  twh_read_labeled_integer(lfailure, FFUNIT,
     &         goodstring, .false.,idname)
         if ( lfailure ) return
      endif
      if ( nangletypes + angle_adj .gt. TBENMAX ) then
         write(6,*) 'RWFORCEFIELD: Number of Angle Terms exceeds '
     &        ,'TBENMAX'
         write(6,*) 'nangletypes,angle_adj TBENMAX',nangletypes
     &        ,angle_adj,TBENMAX
         stop
      endif

      do type = angle_adj+1,angle_adj+nangletypes
         call twh_rw_angle_data(lfailure,rwflag,type,iversion,angle_adj)
         if ( lfailure ) return
      enddo

c     --- input the regular torsion force field information
      goodstring = 'Number of Torsion Terms'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a23)') goodstring
         write(FFUNIT,*) ndihedtypes
      else
         ndihedtypes =  twh_read_labeled_integer(lfailure,FFUNIT,
     &         goodstring, .false.,idname)
         if ( lfailure ) return
      endif
      if ( ndihedtypes + dihed_adj .gt. TTORMAX ) then
         write(6,*) 'RWFORCEFIELD: Number of Torsion Terms exceeds '
     &        ,'TTORMAX'
         write(6,*) 'ndihedtypes,dihed_adj TTORMAX',ndihedtypes
     &        ,dihed_adj,TTORMAX
         stop
      endif

      do type = dihed_adj+1,dihed_adj+ndihedtypes
         call twh_rw_dihedral_data(lfailure,rwflag,type,iversion
     &        ,dihed_adj)
         if ( lfailure ) return
      enddo

c     --- improper torsion force field information
      goodstring = 'Number of Improper Terms'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a24)') goodstring
         write(FFUNIT,*) nimprotypes
      else
         nimprotypes =  twh_read_labeled_integer(lfailure, FFUNIT,
     &         goodstring, .false.,idname)
         if ( lfailure ) return
      endif
      if ( nimprotypes + impr_adj .gt. TIMPMAX ) then
         write(6,*) 'RWFORCEFIELD: Number of Improper Terms exceeds '
     &        ,'TIMPMAX'
         write(6,*) 'nimprotypes,impr_adj TIMPMAX',nimprotypes
     &        ,impr_adj,TIMPMAX
         stop
      endif
      do type = impr_adj+1,impr_adj+nimprotypes
         call twh_rw_improper_data(lfailure, rwflag, type, impr_adj)
         if ( lfailure ) return
      enddo

c     --- angle-angle force field information
      goodstring = 'Number of Angle-Angle Terms'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a27)') goodstring
         write(FFUNIT,*) naatypes
      else
         naatypes =  twh_read_labeled_integer(lfailure, FFUNIT,
     &         goodstring, .false.,idname)
         if ( lfailure ) return
      endif
      if ( naatypes + aa_adj .gt. TAAMAX ) then
         write(6,*) 'RWFORCEFIELD: Number of Angle-Angle Terms exceeds '
     &        ,'TAAMAX'
         write(6,*) 'naatypes,aa_adju TAAMAX',naatypes,aa_adj,TAAMAX
         stop
      endif
      do type = aa_adj+1,aa_adj+naatypes
         call twh_rw_angle_angle_data(lfailure,rwflag, type, aa_adj)
         if ( lfailure ) return
      enddo

c     --- one-five special interactions
      if ( iversion .lt. 13 ) then
c        --- was a mistake in the old text string
         goodstring = 'Number One-Five Types'
      else
         goodstring = 'Number of One-Five Types'
      endif
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a24)') goodstring
         write(FFUNIT,*) noftypes
      else
         if ( iversion .ge. 8 ) then
            noftypes =  twh_read_labeled_integer(lfailure, FFUNIT,
     &         goodstring, .false.,idname)
            if ( lfailure ) return
         else
c           --- set noftypes to 0
            noftypes = 0
         endif
      endif

      if ( noftypes + of_adj .gt. TOFMAX ) then
         write(6,*) 'RWFORCEFIELD: Number of One-Five Types exceeds '
     &        ,' TOFMAX'
         write(6,*) 'noftypes,of_adj TOFMAX',noftypes,of_adj,TOFMAX
         stop
      endif
      do type = of_adj+1,of_adj+noftypes
         call twh_rw_one_five_data(lfailure,rwflag, type, of_adj)
         if ( lfailure ) return
      enddo

c     --- bond increments
      goodstring = 'Number of Bond Increments'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a25)') goodstring
         write(FFUNIT,*) nbitypes
      else
         if ( iversion .ge. 7 ) then
            nbitypes =  twh_read_labeled_integer(lfailure, FFUNIT,
     &           goodstring, .false.,idname)
            if ( lfailure ) return
         else
c           --- set nbitypes to 0
            nbitypes = 0
         endif
      endif
      if ( nbitypes + bi_adj .gt. TVIBMAX ) then
         write(6,*) 'RWFORCEFIELD: Number of Bond Increments exceeds '
     &        ,' TVIBMAX'
         write(6,*) 'nbitypes,bi_adj TVIBMAX',nbitypes,bi_adj,TVIBMAX
         stop
      endif
      do type = bi_adj+1,bi_adj+nbitypes
         call twh_rw_bond_increment_data(lfailure,rwflag,iversion,type
     &        , bi_adj)
         if ( lfailure ) return
      enddo

c     --- H-bond data structures (if needed)
      if ( lrwhbond ) then
         goodstring = 'Number of H-bond Types'
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(a22)') goodstring
         else
            call twh_read_label(lfailure,FFUNIT, goodstring, .false.
     &           ,idname)
            if ( lfailure ) return
         endif
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,*) nhbondtypes
         else
            nhbondtypes = twh_read_integer(lfailure,FFUNIT,
     &            goodstring,.false.)
            if ( lfailure ) return
         endif
         if ( nhbondtypes + hb_adj .gt. THBONDMAX ) then
            write(6,*) 'RWFORCEFIELD: Number of H-bond Types exceeds '
     &           ,' THBONDMAX'
            write(6,*) 'nhbondtypes,hb_adj,THBONDMAX',nhbondtypes,hb_adj
     &           , THBONDMAX
            stop
         endif

         do type = hb_adj + 1, hb_adj + nhbondtypes
            call twh_rw_hbond_data(lfailure,rwflag,type,hb_adj
     &           ,hbcoeffstart,hbcoeffstop)
         enddo
      endif

c     --- close the forcefield file
      rewind(17)
      close(17)

c     --- reset the variables passed in
      ntypes = ntypes + nb_adj
      nbondtypes = nbondtypes + bond_adj
      nbitypes = nbitypes + bi_adj
      nangletypes = nangletypes + angle_adj
      ndihedtypes = ndihedtypes + dihed_adj
      nimprotypes = nimprotypes + impr_adj
      naatypes = naatypes + aa_adj
      noftypes = noftypes + of_adj
      nhbondtypes = nhbondtypes + hb_adj

      return

c     --- error jump out
100   write(6,*) "RWFORCEFIELD: Error opening file ",fffilename
      stop

      end

      subroutine twh_rwtorloop(lfailure,type,rwflag)
c     ******************************************************************
c     * this subroutine reads or writes the ntorloop variable          *
c     *                                                                *
c     * originally written 08-20-2004 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NTORLOOP
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from this subroutine
      logical lfailure
      integer type,rwflag
c     --- local variables
      integer itemp
      character*50 goodstring
c     --- read or write ntorloop
      goodstring = 'Number of Torsion Loops'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a23)') goodstring
         write(FFUNIT,*) twh_get_ntorloop(type)
      else
         itemp = twh_read_labeled_integer(lfailure, FFUNIT,goodstring
     &        , .false.,'rwtorloop')
         if ( lfailure ) return
         call twh_ntorloop(GLB_SET,type,itemp)
      endif
      return
      end


      subroutine setcrossdens(itype)
c     ******************************************************************
c     * this subroutine sets all of the cross densities to the self    *
c     * density.  Useful for reading in old forcefield formats and     *
c     * also for building traditional EAM style forcefields            *
c     *                                                                *
c     * originally written 04-08-2005 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer itype
c     --- local variables
      character*20 eamdensstyle
      integer jtype,itab,icount,eam_ndens
      double precision dvalue

      do jtype = 1,MAXTABTYPE
         if ( itype .ne. jtype ) then
            call twh_eam_ndens(GLB_GET,itype,itype,eam_ndens)
            call twh_eam_ndens(GLB_SET,jtype,itype,eam_ndens)
c           --- set cross terms to the self value
            call twh_eam_dens_style(GLB_GET,itype,itype,eamdensstyle)
            call twh_eam_dens_style(GLB_SET,jtype,itype,eamdensstyle)
            do itab = 1,eam_ndens
               do icount = 1,2
                  call twh_eam_dens(GLB_GET,icount,itype,itype,itab
     &                 ,dvalue)
                  call twh_eam_dens(GLB_SET,icount,jtype,itype,itab
     &                 ,dvalue)
               enddo
            enddo
         endif
      enddo
      return
      end


      subroutine twh_rw_atom_data(lfailure,rwflag, type, nb_adj, ntypes
     &     , iffnum,iversion,lexplicit,lgetcoeffnum,localptype,ncoeff)
c     ******************************************************************
c     * Arguments modified by this subroutine:                         *
c     *  ncoeff                                                        *
c     * Global variables modified by this subroutine::                 *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_BONDPATT
#define FUNCTION_GET_CANAME
#define FUNCTION_GET_CBNAME
#define FUNCTION_GET_CTNAME
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_MASS
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_GET_NBFF
#define FUNCTION_GET_NBNAME
#define FUNCTION_GET_QBASEVALUE
#define FUNCTION_READ_FLOAT
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_SAFE_DOUBLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer rwflag, type, nb_adj, ntypes, iffnum,ncoeff,iversion
      logical lfailure,lexplicit,lgetcoeffnum
      character*30 localptype
c     --- local variables
      logical ltable
c     --- integer scalars
      integer ij,icoeff,jtype,jmax,testtype,itest,jtest,iscan,tabnpair
      character*1 conetemp
      character*2 ename
      character*5 bptemp
      character*(FFNAMELEN) ctemp
      character*15 idname
      character*20 eampairstyle
      character*50 goodstring
c     --- double precision scalars
      double precision dtemp,tpairone,tpairtwo

      idname = 'rwforcefield'
      goodstring = 'Atom Type Number'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring
         write(FFUNIT,*) type
      else
         testtype =  twh_read_labeled_integer(lfailure,FFUNIT
     &        ,goodstring,.false.,idname)
         if ( lfailure ) return
         if ( type-nb_adj .ne. testtype ) then
            write(6,*) 'RWFORCEFIELD: mismatch of Atom types'
            stop
         endif
      endif
      if ( lexplicit ) then
         jmax = ntypes
         if ( iffnum .gt. 1 ) then
            write(6,*) 'RWFORCEFIELD: cannot read multiple ff files'
            write(6,*) 'with explicit nonbonded cross terms'
            stop
         endif
      else
         jmax = type
      endif

      do jtype = type,jmax
         ij = twh_getnbtype(type,jtype)
c           --- initialize ltable
         ltable = .false.
         if ( lgetcoeffnum ) then
c              --- potential with a variable number of nbcoeffs
            if ( localptype .eq. 'Embedded Atom Method' .or.
     &           localptype .eq. 'EAM pair only' ) then
c              --- eam_pair_style
               goodstring = 'eam_pair_style'
               if ( rwflag .eq. WRITE_FLAG ) then
                  write(FFUNIT,'(a14)') goodstring
                  call twh_eam_pair_style(GLB_GET,ij,eampairstyle)
                  write(FFUNIT,'(a20)') eampairstyle
               else
                  call twh_read_labeled_string(lfailure,FFUNIT
     &                 ,goodstring,.false.,idname,eampairstyle)
                  if ( lfailure ) return
                  call twh_eam_pair_style(GLB_SET,ij,eampairstyle)
               endif
c              --- determine number of nbcoeffs (or equiv)
               call twh_eam_pair_style(GLB_GET,ij,eampairstyle)
               if ( eampairstyle .eq. 'table' ) then
                  ltable = .true.
                  ncoeff = 0
               elseif ( eampairstyle .eq. 'exponential' ) then
                  ncoeff = 2
               elseif ( eampairstyle .eq. 'morse' ) then
                  ncoeff = 3
               elseif ( eampairstyle .eq. 'Ackland 3-part' ) then
                  ncoeff = 9
                  ltable = .true.
               elseif ( eampairstyle .eq. 'Ackland Power' ) then
                  ncoeff = 16
                  ltable = .true.
               elseif ( eampairstyle .eq. 'Belashchenko 5-8' ) then
                  ncoeff = 55
               elseif ( eampairstyle .eq. 'none' ) then
                  ncoeff = 0
               else
                  write(6,*) 'rwforcefield: Unknown eam_pair_style'
                  write(6,*) eampairstyle
                  stop
               endif
            else
c                 --- generic tabular potential or multiwell
               ltable = .true.
               ncoeff = 0
            endif
         endif

c           --- nonbonded coefficients or tabular output
         if ( ltable ) then
c           --- table_pair npair
            goodstring = 'table_pair'
            if ( rwflag .eq. WRITE_FLAG ) then
               write(FFUNIT,'(a10)') goodstring
               call twh_table_npair(GLB_GET,ij,tabnpair)
               write(FFUNIT,*) type,jtype,tabnpair
            else
               call twh_read_label(lfailure,FFUNIT,goodstring, .false.
     &              ,idname)
               if ( lfailure ) return
               read(FFUNIT,*) itest,jtest,tabnpair
               call twh_table_npair(GLB_SET,ij,tabnpair)
c                 --- safety check
               if ( itest .ne. type .or. jtest .ne. jtype ) then
                  write(6,*) 
     &                 'RWFORCEFIELD: mismatch of table_pair types'
                  write(6,*) 'type,itest,jtype,jtest',type,itest
     &                 ,jtype,jtest
                  stop
               endif
               if ( tabnpair .gt. MAXTABLE ) then
                  write(6,*) 'RWFORCEFIELD: table_npair > MAXTABLE'
                  write(6,*) 'type,jtype,table_npair,MAXTABLE'
     &                 ,type,jtype,tabnpair,MAXTABLE
                  stop
               endif
            endif

c           --- table_pair data
            goodstring = 'table_pair_data'
            if ( rwflag .eq. WRITE_FLAG ) then
               write(FFUNIT,'(a15)') goodstring
            else
               call twh_read_label(lfailure,FFUNIT,goodstring, .false.
     &              ,idname)
               if ( lfailure ) return
            endif
            call twh_table_npair(GLB_GET,ij,tabnpair)
            do iscan = 1,tabnpair
               if ( rwflag .eq. WRITE_FLAG ) then
                  call twh_table_pair(GLB_GET,1,ij,iscan,tpairone)
                  call twh_table_pair(GLB_GET,2,ij,iscan,tpairtwo)
                  write(FFUNIT,'(e20.10,1x,e20.10)')
     &                 tpairone,tpairtwo
               else
                  read(FFUNIT,*) tpairone,tpairtwo
                  call twh_table_pair(GLB_SET,1,ij,iscan,tpairone)
                  call twh_table_pair(GLB_SET,2,ij,iscan,tpairtwo)
               endif
            enddo
         endif
         if ( ncoeff .gt. 0 ) then
c           --- functional form, coefficients in nbcoeff
            goodstring = 'Nonbond Coefficients'
            if ( rwflag .eq. WRITE_FLAG ) then
               write(FFUNIT,'(a20)') goodstring
            else
               call twh_read_label(lfailure,FFUNIT,goodstring, .false.
     &              ,idname)
               if ( lfailure ) return
            endif
            do icoeff = 1,ncoeff
               if ( rwflag .eq. WRITE_FLAG ) then
                  write(FFUNIT,'(e20.10)') 
     &                 twh_safe_double(twh_get_nbcoeff(icoeff,ij))
               else
                  dtemp = twh_read_float(lfailure,FFUNIT,goodstring
     &                 , .false.)
                  if ( lfailure ) return
                  call twh_nbcoeff(GLB_SET,icoeff,ij,dtemp)
               endif
            enddo
         endif
      enddo

      if ( localptype .eq. 'Embedded Atom Method' .or.
     &     localptype .eq. 'EAM pair only' ) then
         call twh_rw_eam_params(lfailure,rwflag, type, ntypes, iversion)
         if ( lfailure ) return
      endif

c     --- mass
      goodstring = 'Mass'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a4)') goodstring
         write(FFUNIT,'(e20.10)') twh_get_mass(type)
      else
         dtemp = twh_read_labeled_float(lfailure,FFUNIT,goodstring
     &        ,.false.,idname)
         if ( lfailure ) return
         call twh_mass(GLB_SET,type,dtemp)
      endif
c     --- element
      goodstring = 'Element'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a7)') goodstring
         write(FFUNIT,'(a2)') twh_get_elementname(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring,.false.
     &        ,idname,ename)
         if ( lfailure ) return
         call twh_elementname(GLB_SET,type,ename)
      endif
c     --- bond pattern
      goodstring = 'Bond Pattern'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a12)') goodstring
         write(FFUNIT,'(a5)') twh_get_bondpatt(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring
     &        , .false. ,idname,bptemp)
         if ( lfailure ) return
         call twh_bondpatt(GLB_SET,type,bptemp)
      endif
c     --- qbasevalue
      goodstring = 'Base Charge'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a12)') goodstring
         write(FFUNIT,'(e20.10)') twh_get_qbasevalue(type)
      else
         if ( iversion .ge. 7 ) then
            dtemp = twh_read_labeled_float(lfailure, FFUNIT,goodstring
     &           , .false. ,idname)
            if ( lfailure ) return
            call twh_qbasevalue(GLB_SET,type,dtemp)
         else
c              --- older versions did not have this parameter
c              --- just set to zero
            call twh_qbasevalue(GLB_SET,type,0.0d0)
         endif
      endif
c     --- polarizability
      goodstring = 'Polarizability'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a14)') goodstring
         call twh_polarizability(GLB_GET,type,dtemp)
         write(FFUNIT,'(e20.10)') dtemp
      else
         if ( iversion .ge. 13 ) then
            dtemp = twh_read_labeled_float(lfailure,FFUNIT,goodstring
     &           ,.false. ,idname)
            if ( lfailure ) return
            call twh_polarizability(GLB_SET,type,dtemp)
         else
c           --- previous versions did not have polarizability, ignore
         endif
      endif
      if ( localptype .eq. 'Buffered 14-7' ) then
c        --- donor-acceptor flag
         goodstring = 'Donor/Acceptor'
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(a14)') goodstring
            call twh_donor_acceptor(GLB_GET,type,conetemp)
            write(FFUNIT,'(a1)') conetemp
         else
            call twh_read_labeled_string(lfailure,FFUNIT,goodstring
     &           ,.false. ,idname,conetemp)
            if ( lfailure ) return
            call twh_donor_acceptor(GLB_SET,type,conetemp)
         endif
      endif
c     --- force field name
      goodstring = 'Force Field Name'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring
         write(FFUNIT,'(a10)') twh_get_nbff(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring,.false.
     &        ,idname,ctemp)
         if ( lfailure ) return
         call twh_nbff(GLB_SET,type,ctemp)
      endif
c        --- input the names of those atoms
      goodstring = 'Atom Names'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring, .false.,idname)
         if ( lfailure ) return
      endif
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') twh_get_nbname(type)
         write(FFUNIT,'(a10)') twh_get_cbname(type)
         write(FFUNIT,'(a10)') twh_get_caname(type)
         write(FFUNIT,'(a10)') twh_get_ctname(type)
      else
         call twh_read_string(lfailure,FFUNIT,goodstring,.false.,ctemp)
         if ( lfailure ) return
         call twh_nbname(GLB_SET,type,ctemp)
         call twh_read_string(lfailure,FFUNIT,goodstring,.false.,ctemp)
         if ( lfailure ) return
         call twh_cbname(GLB_SET,type,ctemp)
         call twh_read_string(lfailure,FFUNIT,goodstring,.false.,ctemp)
         if ( lfailure ) return
         call twh_caname(GLB_SET,type,ctemp)
         call twh_read_string(lfailure,FFUNIT,goodstring,.false.,ctemp)
         if ( lfailure ) return
         call twh_ctname(GLB_SET,type,ctemp)
      endif

      return
      end


      subroutine twh_rw_bond_data(lfailure,iversion,rwflag,type
     &     , bond_adj)
c     ******************************************************************
c     * Arguments modified by this subroutine:                         *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 07-14-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NVSAME
#define FUNCTION_GET_VIBCOEFF
#define FUNCTION_GET_VIBFF
#define FUNCTION_GET_VIBNAMES
#define FUNCTION_READ_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_SAFE_DOUBLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer iversion, rwflag, type, bond_adj
c     --- local variables
c     --- integer scalars
      integer coeffstart,coeffstop,isame,icoeff,iatom,vstemp,testtype
c     --- character string scalars
      character*(FFNAMELEN) order
      character*15 idname
      character*50 goodstring
c     --- character string arrays
      character*(FFNAMELEN) ctemp
      dimension ctemp(MAXVNAME)
c     --- double precision scalars
      double precision dtemp

      idname = 'rwforcefield'
      goodstring = 'Bond Type Number'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring
         write(FFUNIT,*) type
      else
         testtype =  twh_read_labeled_integer(lfailure, FFUNIT,
     &      goodstring, .false.,idname)
         if ( lfailure ) return
         if ( type-bond_adj .ne. testtype ) then
            write(6,*) 'RWFORCEFIELD: mismatch of Bond types'
            stop
         endif
      endif

      goodstring = 'Bond Style'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') goodstring
         call twh_vibstyle(GLB_GET,type,vstemp)
         write(FFUNIT,*) vstemp
      else
         vstemp =  twh_read_labeled_integer(lfailure, FFUNIT,
     &        goodstring, .false.,idname)
         if ( lfailure ) return
         call twh_vibstyle(GLB_SET,type,vstemp)
      endif

      call twh_get_vibration_params(iversion,type,vstemp,coeffstart
     &     , coeffstop)

      goodstring = 'Bond Coefficients'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a17)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring, .false.,idname)
         if ( lfailure ) return
      endif
      do icoeff = coeffstart,coeffstop
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(e20.10)') 
     &           twh_safe_double(twh_get_vibcoeff(type,icoeff))
         else
            dtemp = twh_read_float(lfailure,FFUNIT,goodstring,.false.)
            if ( lfailure ) return
            call twh_vibcoeff(GLB_SET,type,icoeff,dtemp)
         endif
      enddo

c     --- vibration order for this type of interaction
      goodstring = 'Vibration Order'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a15)') goodstring
         call twh_viborder(GLB_GET,type,order)
         write(FFUNIT,'(a10)') order
      else
         if ( iversion .ge. 14 ) then
c           --- read in the viborder
            call twh_read_label(lfailure,FFUNIT,goodstring, .false.
     &           ,idname)
            if ( lfailure ) return
            read(FFUNIT,'(a10)') order
         else
c           --- prior versions did not have viborder; set to wild
            order = 'wild'
         endif
         call twh_viborder(GLB_SET,type,order)
      endif

c     --- force field name
      goodstring = 'Force Field Name'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring
         write(FFUNIT,'(a10)') twh_get_vibff(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring
     &        , .false.,idname,ctemp(1))
         if ( lfailure ) return
         call twh_vibff(GLB_SET,type,ctemp(1))
      endif

c     --- number of atoms in this force field that
c     --- use the same parameters
      goodstring = 'Number of Atoms with Same Parameters'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a38)') goodstring
         write(FFUNIT,*) twh_get_nvsame(type)
      else
         isame =  twh_read_labeled_integer(lfailure, FFUNIT,
     &      goodstring, .false.,idname)
         if ( lfailure ) return
         call twh_nvsame(GLB_SET,type,isame)
      endif
c     --- names of those atoms
      goodstring = 'Atom Names'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      do isame = 1,twh_get_nvsame(type)
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(2(a10,1x))') 
     &           (twh_get_vibnames(type,isame,iatom),iatom=1,MAXVNAME)
         else
            read(FFUNIT,'(2(a10,1x))') (ctemp(iatom),iatom=1,MAXVNAME)
            do iatom = 1,MAXVNAME
               call twh_vibnames(GLB_SET,type,isame,iatom,ctemp(iatom))
            enddo
         endif
      enddo

      return
      end


      subroutine twh_rw_angle_data(lfailure,rwflag,type,iversion
     &     ,angle_adj)
c     ******************************************************************
c     * Arguments modified by this subroutine:                         *
c     * Global variables modified by this subroutine::                 *
c     *  anglestyle                                                    *
c     *  bendnames                                                     *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ANGLESTYLE
#define FUNCTION_GET_BENCOEFF
#define FUNCTION_GET_BENDFF
#define FUNCTION_GET_BENDNAMES
#define FUNCTION_GET_LBONANG
#define FUNCTION_GET_LBONBON
#define FUNCTION_GET_NBSAME
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_READ_LABELED_LOGICAL
#define FUNCTION_READ_FLOAT
#define FUNCTION_SAFE_DOUBLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag, type, iversion
      integer angle_adj
c     --- local variables
      logical lrwbondangle,lrwbondbond,lbastring,lreadtemp
      integer coeffstart,coeffstop
      integer isame,icoeff,iatom,itemp
      integer testtype
      integer ibastart,ibaend,ibbstart,ibbend
      character*(FFNAMELEN) ctemp,ctemparray
      dimension ctemparray(MAXBNAME)
      character*15 idname,order
      character*50 goodstring
      double precision dtemp

      idname = 'rwforcefield'
      goodstring = 'Angle Type Number'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a17)') goodstring
         write(FFUNIT,*) type
      else
         testtype =  twh_read_labeled_integer(lfailure, FFUNIT,
     &      goodstring, .false.,idname)
         if ( lfailure ) return
         if ( type-angle_adj .ne. testtype ) then
            write(6,*) 'RWFORCEFIELD: mismatch of Angle types'
            stop
         endif
      endif

c     --- initialize lrwbondangle and lrwbondbond to most common 
c     --- value of false.  Then we can set these to true below if 
c     --- needed
      lrwbondangle = .false.
      lrwbondbond = .false.

      goodstring = 'Angle Style'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a11)') goodstring
         write(FFUNIT,*) twh_get_anglestyle(type)
      else
         itemp =  twh_read_labeled_integer(lfailure, FFUNIT,goodstring
     &        , .false.,idname)
         if ( lfailure ) return
         call twh_anglestyle(GLB_SET,type,itemp)
      endif

      itemp = twh_get_anglestyle(type)
      call get_angle_params(itemp, iversion,coeffstart, coeffstop, 
     &   ibastart, ibaend, ibbstart, ibbend,
     &   lrwbondangle, lbastring, lrwbondbond)

c        --- output bond-angle if needed
      if ( lrwbondangle ) then
         goodstring = 'Bond-Angle Logical'
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(a18)') goodstring
            write(FFUNIT,*) twh_get_lbonang(type)
         else
            lreadtemp =  twh_read_labeled_logical(lfailure, FFUNIT,
     &         goodstring, .false.,idname)
            if ( lfailure ) return
            call twh_lbonang(GLB_SET,type,lreadtemp)
         endif
         if ( twh_get_lbonang(type) ) then
            if ( iversion .lt. 14 ) then
               if ( lbastring ) then
                  if ( rwflag .eq. READ_FLAG ) then
c                    --- skip three lines as this is no longer needed
                     read(FFUNIT,*)
                     read(FFUNIT,*)
                     read(FFUNIT,*)
                  endif
               endif
            endif
            goodstring = 'Bond-Angle Coeffecients'
            if ( rwflag .eq. WRITE_FLAG ) then
               write(FFUNIT,'(a23)') goodstring
            else
               call twh_read_label(lfailure,FFUNIT,goodstring,.false.
     &              ,idname)
               if ( lfailure ) return
            endif
            do icoeff = ibastart,ibaend
               if ( rwflag .eq. WRITE_FLAG ) then
                  write(FFUNIT,'(e20.10)') 
     &                 twh_safe_double(twh_get_bencoeff(type,icoeff))
               else
                  dtemp = twh_read_float(lfailure,FFUNIT, goodstring
     &                 ,.false.)
                  if ( lfailure ) return
                  call twh_bencoeff(GLB_SET,type,icoeff,dtemp)
               endif
            enddo
         endif
      endif

c        --- output bond-bond if needed
      if ( lrwbondbond ) then
         goodstring = 'Bond-Bond Logical'
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(a17)') goodstring
            write(FFUNIT,*) twh_get_lbonbon(type)
         else
            lreadtemp =  twh_read_labeled_logical(lfailure, FFUNIT,
     &            goodstring, .false.,idname)
            if ( lfailure ) return
            call twh_lbonbon(GLB_SET,type,lreadtemp)
         endif
         if ( twh_get_lbonbon(type) ) then
            goodstring = 'Bond-Bond Coeffecients'
            if ( rwflag .eq. WRITE_FLAG ) then
               write(FFUNIT,'(a22)') goodstring
            else
               call twh_read_label(lfailure,FFUNIT,goodstring,.false.
     &              ,idname)
               if ( lfailure ) return
            endif
            do icoeff = ibbstart,ibbend
               if ( rwflag .eq. WRITE_FLAG ) then
                  write(FFUNIT,'(e20.10)') twh_get_bencoeff(type,icoeff)
               else
                  dtemp = twh_read_float(lfailure,FFUNIT, goodstring
     &                 ,.false.)
                  if ( lfailure ) return
                  call twh_bencoeff(GLB_SET,type,icoeff,dtemp)
               endif
            enddo
         endif
      endif

      goodstring = 'Angle Coefficients'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a18)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring, .false.,idname)
         if ( lfailure ) return
      endif
      do icoeff = coeffstart,coeffstop
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(e20.10)') 
     &           twh_safe_double(twh_get_bencoeff(type,icoeff))
         else
            dtemp = twh_read_float(lfailure,FFUNIT, goodstring,.false.)
            if ( lfailure ) return
            call twh_bencoeff(GLB_SET,type,icoeff,dtemp)
         endif
      enddo

c     --- angle order for this type of interaction
      goodstring = 'Angle Order'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a11)') goodstring
         call twh_angleorder(GLB_GET,type,order)
         write(FFUNIT,'(a15)') order
      else
         if ( iversion .ge. 14 ) then
c           --- read in the angleorder
            call twh_read_label(lfailure,FFUNIT,goodstring,.false.
     &           ,idname)
            if ( lfailure ) return
            read(FFUNIT,'(a15)') order
         else
c           --- prior versions did not have angleorder; set to wild
            order = 'wild'
         endif
         call twh_angleorder(GLB_SET,type,order)
      endif

c     --- force field name
      goodstring = 'Force Field Name'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring
         write(FFUNIT,'(a10)') twh_get_bendff(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring
     &        , .false.,idname,ctemp)
         if ( lfailure ) return
         call twh_bendff(GLB_SET,type,ctemp)
      endif

c        --- number of atoms in this force field that
c        --- use the same parameters
      goodstring = 'Number of Atoms with Same Parameters'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a38)') goodstring
         write(FFUNIT,*) twh_get_nbsame(type)
      else
         itemp = twh_read_labeled_integer(lfailure,FFUNIT,goodstring
     &        ,.false.,idname)
         if ( lfailure ) return
         call twh_nbsame(GLB_SET,type,itemp)
      endif
c        --- names of those atoms
      goodstring = 'Atom Names'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT, goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      do isame = 1,twh_get_nbsame(type)
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(3(a10,1x))') 
     &           (twh_get_bendnames(type,isame,iatom),iatom=1,MAXBNAME)
         else
            read(FFUNIT,'(3(a10,1x))') 
     &           (ctemparray(iatom),iatom=1,MAXBNAME)
            do iatom = 1,MAXBNAME
               call twh_bendnames(GLB_SET,type,isame,iatom
     &              ,ctemparray(iatom))
            enddo
         endif
      enddo

      return
      end


      subroutine twh_rw_dihedral_data(lfailure,rwflag,type,iversion
     &     , dihed_adj)
c     ******************************************************************
c     * Arguments modified by this subroutine:                         *
c     * Global variables modified by this subroutine::                 *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LOFTOR
#define FUNCTION_GET_NTSAME
#define FUNCTION_GET_TORCOEFF
#define FUNCTION_GET_TORFF
#define FUNCTION_GET_TORNAMES
#define FUNCTION_GET_TORSTRING
#define FUNCTION_GET_TORSTYLE
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_LOGICAL
#define FUNCTION_READ_FLOAT
#define FUNCTION_SAFE_DOUBLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag, type, iversion, dihed_adj
c     --- local variables
      logical ltorstring,lreadtemp
      integer isame,icoeff,iatom,itemp
      integer coeffstart,coeffstop
      integer testtype
      character*(FFNAMELEN) ctemp
      dimension ctemp(MAXTNAME)
      character*15 idname,order
      character*50 goodstring
      double precision dtemp

      idname = 'rwforcefield'
      goodstring = 'Torsion Type Number'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a19)') goodstring
         write(FFUNIT,*) type
      else
         testtype =  twh_read_labeled_integer(lfailure,FFUNIT,
     &         goodstring, .false.,idname)
         if ( lfailure ) return
         if ( type-dihed_adj .ne. testtype ) then
            write(6,*) 'RWFORCEFIELD: mismatch of Torsion types'
            stop
         endif
      endif

      goodstring = 'Torsion Style'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a13)') goodstring
         write(FFUNIT,*) twh_get_torstyle(type)
      else
         itemp = twh_read_labeled_integer(lfailure,FFUNIT,
     &      goodstring, .false.,idname)
         if ( lfailure ) return
         call twh_torstyle(GLB_SET,type,itemp)
      endif

      goodstring = 'One-Four Nonbond Logical'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a24)') goodstring
         write(FFUNIT,*) twh_get_loftor(type)
      else
         lreadtemp =  twh_read_labeled_logical(lfailure, FFUNIT,
     &      goodstring, .false.,idname)
         if ( lfailure ) return
         call twh_loftor(GLB_SET,type,lreadtemp)
      endif
      if ( twh_get_loftor(type) ) then
         goodstring = 'One-Four Coulombic Scaling'
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(a26)') goodstring
            write(FFUNIT,'(e20.10)') 
     &           twh_safe_double(twh_get_torcoeff(type,-1))
         else
            dtemp = twh_read_labeled_float(lfailure,FFUNIT,goodstring
     &           , .false. ,idname)
            if ( lfailure ) return
            call twh_torcoeff(GLB_SET,type,-1,dtemp)
         endif
      endif

      call twh_get_torsion_params(lfailure,iversion,rwflag,type
     &     ,coeffstart, coeffstop, ltorstring)
      if ( lfailure ) return

      if ( ltorstring ) then
         goodstring = 'Torsion String Names'
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(a20)') goodstring
         else
            call twh_read_label(lfailure,FFUNIT, goodstring, .false.
     &           ,idname)
            if ( lfailure ) return
         endif
         do icoeff = 1,4
            if ( rwflag .eq. WRITE_FLAG ) then
               write(FFUNIT,'(a10,1x)') twh_get_torstring(type,icoeff)
            else
               call twh_read_string(lfailure,FFUNIT,goodstring,.false.
     &              ,ctemp(1))
               if ( lfailure ) return
               call twh_torstring(GLB_SET,type,icoeff,ctemp(1))
            endif
         enddo
      endif
      goodstring = 'Torsion Coefficients'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a20)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      do icoeff = coeffstart,coeffstop
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(e20.10)') 
     &           twh_safe_double(twh_get_torcoeff(type,icoeff))
         else
            dtemp = twh_read_float(lfailure,FFUNIT,goodstring,.false.)
            if ( lfailure ) return
            call twh_torcoeff(GLB_SET,type,icoeff,dtemp)
         endif
      enddo

c     --- torsion order for this type of interaction
      goodstring = 'Torsion Order'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a13)') goodstring
         call twh_torsionorder(GLB_GET,type,order)
         write(FFUNIT,'(a15)') order
      else
         if ( iversion .ge. 14 ) then
c           --- read in the torsionorder
            call twh_read_label(lfailure,FFUNIT,goodstring,.false.
     &           ,idname)
            if ( lfailure ) return
            read(FFUNIT,'(a15)') order
         else
c           --- prior versions did not have torsionorder; set to wild
            order = 'wild'
         endif
         call twh_torsionorder(GLB_SET,type,order)
      endif

c     --- force field name
      goodstring = 'Force Field Name'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring
         write(FFUNIT,'(a10)') twh_get_torff(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring
     &        , .false.,idname,ctemp(1))
         if ( lfailure ) return
         call twh_torff(GLB_SET,type,ctemp(1))
      endif
c        --- number of atoms in this force field that
c        --- use the same parameters
      goodstring = 'Number of Atoms with Same Parameters'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a38)') goodstring
         write(FFUNIT,*) twh_get_ntsame(type)
      else
         itemp = twh_read_labeled_integer(lfailure, FFUNIT,goodstring
     &        , .false.,idname)
         if ( lfailure ) return
         call twh_ntsame(GLB_SET,type,itemp)
      endif
c        --- safety check value of ntsame
      if ( twh_get_ntsame(type) .gt. MAXSAMET ) then
c           --- out of bounds
         write(6,*) 'RWFORCEFIELD: ntsame exceeds MAXSAMET'
         write(6,*) 'MAXSAMET (set in preproc.h):',MAXSAMET
         write(6,*) 'type,ntsame',type,twh_get_ntsame(type)
         stop
      endif
c        --- names of those atoms
      goodstring = 'Atom Names'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      do isame = 1,twh_get_ntsame(type)
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,22)
     &           (twh_get_tornames(type,isame,iatom),iatom=1,MAXTNAME)
         else
            read(FFUNIT,22)
     &           (ctemp(iatom),iatom=1,MAXTNAME)
            do iatom = 1,MAXTNAME
               call twh_tornames(GLB_SET,type,isame,iatom,ctemp(iatom))
            enddo
         endif
 22      format(MAXTNAME(a10,1x)) 
      enddo

      return
      end


      subroutine twh_rw_improper_data(lfailure, rwflag, type, impr_adj)
c     ******************************************************************
c     * Arguments modified by this subroutine:                         *
c     * Global variables modified by this subroutine::                 *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPCOEFF
#define FUNCTION_GET_IMPFF
#define FUNCTION_GET_IMPFORM
#define FUNCTION_GET_IMPNAMES
#define FUNCTION_GET_IMPSTYLE
#define FUNCTION_GET_NIMPSAME
#define FUNCTION_READ_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag, type, impr_adj
c     --- local variables
      integer isame,icoeff,iatom,itemp
      integer coeffstart,coeffstop
      integer testtype
      character*(FFNAMELEN) ctemp
      dimension ctemp(4)
      character*15 idname
      character*50 goodstring
      double precision dtemp

      idname = 'rwforcefield'
      goodstring = 'Improper Type Number'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a20)') goodstring
         write(FFUNIT,*) type
      else
         testtype =  twh_read_labeled_integer(lfailure, FFUNIT,
     &      goodstring, .false.,idname)
         if ( lfailure ) return
         if ( type-impr_adj .ne. testtype ) then
            write(6,*) 'RWFORCEFIELD: mismatch of Improper types'
            stop
         endif
      endif
      goodstring = 'Improper Form'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a13)') goodstring
         write(FFUNIT,*) twh_get_impform(type)
      else
         itemp = twh_read_labeled_integer(lfailure,FFUNIT,goodstring
     &        ,.false.,idname)
         if ( lfailure ) return
         call twh_impform(GLB_SET,type,itemp)
      endif
      goodstring = 'Improper Style'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a14)') goodstring
         write(FFUNIT,*) twh_get_impstyle(type)
      else
         itemp = twh_read_labeled_integer(lfailure,FFUNIT,goodstring
     &        ,.false.,idname)
         if ( lfailure ) return
         call twh_impstyle(GLB_SET,type,itemp)
      endif

      itemp = twh_get_impstyle(type)
      call get_improper_params(itemp,coeffstart,coeffstop)

      goodstring = 'Improper Coefficients'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a21)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      do icoeff = coeffstart,coeffstop
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(e20.10)') twh_get_impcoeff(type,icoeff)
         else
            dtemp = twh_read_float(lfailure,FFUNIT,goodstring,.false.)
            if ( lfailure ) return
            call twh_impcoeff(GLB_SET,type,icoeff,dtemp)
         endif
      enddo
c        --- force field name
      goodstring = 'Force Field Name'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring 
         write(FFUNIT,'(a10)') twh_get_impff(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring,.false.
     &        ,idname,ctemp(1))
         if ( lfailure ) return
         call twh_impff(GLB_SET,type,ctemp(1))
      endif
c        --- number of atoms in this force field that
c        --- use the same parameters
      goodstring = 'Number of Atoms with Same Parameters'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a38)') goodstring
         write(FFUNIT,*) twh_get_nimpsame(type)
      else
         itemp = twh_read_labeled_integer(lfailure,FFUNIT,goodstring
     &        ,.false.,idname)
         if ( lfailure ) return
         call twh_nimpsame(GLB_SET,type,itemp)
      endif
c        --- names of those atoms
      goodstring = 'Atom Names'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      do isame = 1,twh_get_nimpsame(type)
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(4(a10,1x))') 
     &           (twh_get_impnames(type,isame,iatom),iatom=1,4)
         else
            read(FFUNIT,'(4(a10,1x))') (ctemp(iatom),iatom=1,4)
            do iatom = 1,4
               call twh_impnames(GLB_SET,type,isame,iatom,ctemp(iatom))
            enddo
         endif
      enddo

      return
      end

      subroutine twh_rw_angle_angle_data(lfailure,rwflag, type, aa_adj)
c     ******************************************************************
c     * Arguments modified by this subroutine:                         *
c     * Global variables modified by this subroutine::                 *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_AACOEFF
#define FUNCTION_GET_AAFF
#define FUNCTION_GET_AANAMES
#define FUNCTION_GET_AASTYLE
#define FUNCTION_GET_NAASAME
#define FUNCTION_READ_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from this subroutine
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag, type, aa_adj
c     --- local variables
      integer isame,icoeff,iatom
      integer coeffstart,coeffstop
      integer testtype,itemp
      character*(FFNAMELEN) ctemp
      dimension ctemp(4)
      character*15 idname
      character*50 goodstring
      double precision dtemp

      idname = 'rwforcefield'
      goodstring = 'Angle-Angle Type Number'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a23)') goodstring
         write(FFUNIT,*) type
      else
         testtype =  twh_read_labeled_integer(lfailure, FFUNIT,
     &      goodstring, .false.,idname)
         if ( lfailure ) return
         if ( type-aa_adj .ne. testtype ) then
            write(6,*) 'RWFORCEFIELD: mismatch of Angle-Angle types'
            stop
         endif
      endif
      goodstring = 'Angle-Angle Style'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a17)') goodstring
         write(FFUNIT,*) twh_get_aastyle(type)
      else
         itemp = twh_read_labeled_integer(lfailure, FFUNIT,goodstring
     &        ,.false.,idname)
         if ( lfailure ) return
         call twh_aastyle(GLB_SET,type,itemp)
      endif

      if ( twh_get_aastyle(type) .eq. 1 ) then
c           --- Compass with build in angle terms
         coeffstart = 0
         coeffstop = 0
      elseif ( twh_get_aastyle(type) .eq. 2 ) then
c           --- Compass with explicitly specified angle terms
         coeffstart = 0
         coeffstop = 2
      else
         write(6,*) 'RWFORCEFIELD: invalid aastyle'
     &        ,twh_get_aastyle(type)
         stop
      endif

      goodstring = 'Angle-Angle Coefficients'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a24)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      do icoeff = coeffstart,coeffstop
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(e20.10)') twh_get_aacoeff(type,icoeff)
         else
            dtemp = twh_read_float(lfailure,FFUNIT,goodstring,.false.)
            if ( lfailure ) return
            call twh_aacoeff(GLB_SET,type,icoeff,dtemp)
         endif
      enddo
c        --- force field name
      goodstring = 'Force Field Name'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring
         write(FFUNIT,'(a10)') twh_get_aaff(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring,.false.
     &        ,idname,ctemp(1))
         if ( lfailure ) return
         call twh_aaff(GLB_SET,type,ctemp(1))
      endif
c        --- number of atoms in this force field that
c        --- use the same parameters
      goodstring = 'Number of Atoms with Same Parameters'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a38)') goodstring
         write(FFUNIT,*) twh_get_naasame(type)
      else
         itemp = twh_read_labeled_integer(lfailure,FFUNIT,goodstring
     &        ,.false.,idname)
         if ( lfailure ) return
         call twh_naasame(GLB_SET,type,itemp)
      endif
c        --- names of those atoms
      goodstring = 'Atom Names'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      do isame = 1,twh_get_naasame(type)
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(4(a10,1x))') 
     &           (twh_get_aanames(type,isame,iatom),iatom=1,4)
         else
            read(FFUNIT,'(4(a10,1x))') (ctemp(iatom),iatom=1,4)
            do iatom = 1,4
               call twh_aanames(GLB_SET,type,isame,iatom,ctemp(iatom))
            enddo
         endif
      enddo

      return
      end


      subroutine twh_rw_one_five_data(lfailure, rwflag, type, of_adj)
c     ******************************************************************
c     * Arguments modified by this subroutine:                         *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_OFCOEFF
#define FUNCTION_GET_OFFF
#define FUNCTION_GET_OFNAMES
#define FUNCTION_GET_OFSTYLE
#define FUNCTION_READ_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag,type,itemp,of_adj
c     --- local variables
      integer coeffstart,coeffstop
      integer icoeff,iatom,testtype
      character*(FFNAMELEN) ctemp
      dimension ctemp(MAXOFNAME)
      character*15 idname
      character*50 goodstring
      double precision dtemp

      idname = 'rwforcefield'
      goodstring = 'One-Five Type Number'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a20)') goodstring
         write(FFUNIT,*) type
      else
         testtype =  twh_read_labeled_integer(lfailure, FFUNIT,
     &      goodstring, .false.,idname)
         if ( lfailure ) return
         if ( type - of_adj .ne. testtype ) then
            write(6,*) 'RWFORCEFIELD: mismatch One-Five types'
            stop
         endif
      endif

      goodstring = 'One-Five Style'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a14)') goodstring
         write(FFUNIT,*) twh_get_ofstyle(type)
      else
         itemp = twh_read_labeled_integer(lfailure,FFUNIT,goodstring
     &        ,.false.,idname)
         if ( lfailure ) return
         call twh_ofstyle(GLB_SET,type,itemp)
      endif
      if ( twh_get_ofstyle(type) .eq. 1 ) then
c           --- lennard-jones with special parameters
         coeffstart = 1
         coeffstop = 2
      elseif ( twh_get_ofstyle(type) .eq. 2 ) then
c           --- a/r^12
         coeffstart = 1
         coeffstop = 1
      else
         write(6,*) 'RWFORCEFIELD: invalid ofstyle'
     &        ,twh_get_ofstyle(type)
         stop
      endif

      goodstring = 'One-Five Coefficients'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a22)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      do icoeff = coeffstart,coeffstop
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(e20.10)') twh_get_ofcoeff(type,icoeff)
         else
            dtemp = twh_read_float(lfailure,FFUNIT,goodstring,.false.)
            if ( lfailure ) return
            call twh_ofcoeff(GLB_SET,type,icoeff,dtemp)
         endif
      enddo

c        --- force field name
      goodstring = 'Force Field Name'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring
         write(FFUNIT,'(a10)') twh_get_offf(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring,.false.
     &        ,idname,ctemp(1))
         if ( lfailure ) return
         call twh_offf(GLB_SET,type,ctemp(1))
      endif

c        --- names of those atoms
      goodstring = 'Atom Names'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(5(a10,1x))') 
     &        (twh_get_ofnames(type,iatom),iatom=1,MAXOFNAME)
      else
         read(FFUNIT,'(5(a10,1x))') (ctemp(iatom),iatom=1,MAXOFNAME)
         do iatom = 1,MAXOFNAME
            call twh_ofnames(GLB_SET,type,iatom,ctemp(iatom))
         enddo
      endif

      return
      end


      subroutine twh_rw_bond_increment_data(lfailure,rwflag, iversion
     &     , type, bi_adj)
c     ******************************************************************
c     * Arguments modified by this subroutine:                         *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_QBIFF
#define FUNCTION_GET_QBINAMES
#define FUNCTION_GET_QBIVALUE
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from this subroutine
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag, iversion, type, bi_adj
c     --- local variables
      integer iatom
      integer testtype
      character*FFNAMELEN order
      character*(FFNAMELEN) ctemp
      dimension ctemp(2)
      character*15 idname
      character*50 goodstring
      double precision dtemp

      idname = 'rwforcefield'
      goodstring = 'Bond Increment Type Number'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a26)') goodstring
         write(FFUNIT,*) type
      else
         testtype =  twh_read_labeled_integer(lfailure, FFUNIT,
     &      goodstring, .false.,idname)
         if ( lfailure ) return
         if ( type-bi_adj .ne. testtype ) then
            write(6,*) 'RWFORCEFIELD: mismatch Bond Increment types'
            stop
         endif
      endif

      goodstring = 'Bond Increment Value'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a20)') goodstring
         write(FFUNIT,'(e20.10)') twh_get_qbivalue(type)
      else
         dtemp = twh_read_labeled_float(lfailure,FFUNIT,goodstring
     &        ,.false. ,idname)
         if ( lfailure ) return
         call twh_qbivalue(GLB_SET,type,dtemp)
      endif
c     --- bond increment order
      goodstring = 'Bond Increment Order'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a20)') goodstring
         call twh_biorder(GLB_GET,type,order)
         write(FFUNIT,'(a10)') order
      else
         if ( iversion .ge. 14 ) then
            call twh_read_labeled_string(lfailure,FFUNIT,goodstring
     &           ,.false.,idname,order)
            if ( lfailure ) return
         else
c           --- older versions did not have bond increment order
            order = 'wild'
         endif
         call twh_biorder(GLB_SET,type,order)
      endif

c     --- force field name
      goodstring = 'Force Field Name'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a16)') goodstring
         write(FFUNIT,'(a10)') twh_get_qbiff(type)
      else
         call twh_read_labeled_string(lfailure,FFUNIT,goodstring,.false.
     &        ,idname,ctemp(1))
         if ( lfailure ) return
         call twh_qbiff(GLB_SET,type,ctemp(1))
      endif
c     --- names of those atoms
      goodstring = 'Atom Names'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a10)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(2(a10,1x))') 
     &        (twh_get_qbinames(type,iatom),iatom=1,2)
      else
         read(FFUNIT,'(2(a10,1x))') (ctemp(iatom),iatom=1,2)
         call twh_set_qbinames(type,ctemp)
      endif

      return
      end

      subroutine twh_rw_hbond_data( lfailure,rwflag, type, hb_adj
     &     , hbcoeffstart, hbcoeffstop)
c     ******************************************************************
c     * Arguments modified by this subroutine:                         *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_HBONDCOEFF
#define FUNCTION_GET_HBONDNAMES
#define FUNCTION_READ_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag, type, hb_adj, hbcoeffstart, hbcoeffstop
c     --- local variables
      integer icoeff,iatom
      integer testtype
      character*(FFNAMELEN) ctemp
      dimension ctemp(MAXHBONDNAME)
      character*15 idname
      character*50 goodstring
      double precision dtemp

      idname = 'rwforcefield'
      goodstring = 'H-bond Type Number'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a18)') goodstring
         write(FFUNIT,*) type
      else
         testtype =  twh_read_labeled_integer(lfailure, FFUNIT,
     &       goodstring, .false.,idname)
         if ( lfailure ) return
         if ( type - hb_adj .ne. testtype ) then
            write(6,*) 'RWFORCEFIELD: mismatch H-bond types'
            stop
         endif
      endif

c           --- right now the h-bond style depends upon the 
c           --- classical_potential so there is no reason to list the
c           --- style here.  The values of hbcoeffstart and hbcoeffstop
c           --- are set back up in the classical_potential section
      goodstring = 'H-bond Coefficients'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a19)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif

      do icoeff = hbcoeffstart,hbcoeffstop
         if ( rwflag .eq. WRITE_FLAG ) then
            write(FFUNIT,'(e20.10)') twh_get_hbondcoeff(type,icoeff)
         else
            dtemp = twh_read_float(lfailure,FFUNIT,goodstring,.false.)
            if ( lfailure ) return
            call twh_hbondcoeff(GLB_SET,type,icoeff,dtemp)
         endif
      enddo

      goodstring = 'H-bond Names'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a12)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring,.false.,idname)
         if ( lfailure ) return
      endif
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(2(a10,1x))') 
     &         (twh_get_hbondnames(type,iatom),iatom=1,MAXHBONDNAME)
      else
         read(FFUNIT,'(2(a10,1x))') (ctemp(iatom),iatom=1,MAXHBONDNAME)
         do iatom = 1,MAXHBONDNAME
            call twh_hbondnames(GLB_SET,type,iatom,ctemp(iatom))
         enddo
      endif

      return
      end

      subroutine twh_get_potential_type(lfailure,iversion, localptype)
c     ******************************************************************
c     * Obtain the potential type                                      *
c     * Arguments set in this procedure:                               *
c     *  localptype                                                    *
c     *                                                                *
c     * Error conditions:                                              *
c     *  This function will halt upon an error condition               *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#define FUNCTION_READ_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer iversion
      character*30 localptype
c     --- local variables
      logical lproblem
      integer iptype
      character*15 idname
      character*30 truecp
      character*50 goodstring

      idname = 'rwforcefield'
      goodstring = 'Potential Type'
      call twh_read_label(lfailure,FFUNIT,goodstring, .false.,idname)
      if ( lfailure ) return
      if ( iversion .lt. 6 ) then
c           --- read in old integer code and translate into current
c           --- character string
         iptype = twh_read_integer(lfailure,FFUNIT,goodstring,.false.)
         if ( lfailure ) return
         if ( iptype .eq. -3 ) then
            localptype = 'Repulsive Sphere'
         elseif ( iptype .eq. 0 ) then
            localptype = 'Lennard-Jones'
         elseif ( iptype .eq. 1 ) then
            localptype = '9-6'
         elseif ( iptype .eq. 2 ) then
            localptype = 'Exponential-6'
         elseif ( iptype .eq. 3 ) then
            localptype = 'Hard Sphere'
         elseif ( iptype .eq. 4 ) then
            localptype = 'Exponential-12-6'
         elseif ( iptype .eq. 5 ) then
            localptype = 'Stillinger-Weber'
         elseif ( iptype .eq. 6 ) then
            localptype = 'Embedded Atom Method'
         elseif ( iptype .eq. 7 ) then
            localptype = '12-6 plus solvation'
         elseif ( iptype .eq. 8 ) then
            localptype = '12-9-6'
         elseif ( iptype .eq. 9 ) then
            localptype = 'Square Well'
         elseif ( iptype .eq. 10 ) then
            localptype = 'Tabulated Pair'
         else
            write(6,*) 'RWFORCEFIELD: unknown iptype',iptype
            stop
         endif
      else
c           --- read in the text string
         call twh_read_string(lfailure,FFUNIT,goodstring,.false.
     &        ,localptype)
         if ( lfailure ) return
      endif
      if ( localptype .ne. twh_get_classical_potential() ) then
         lproblem = .true.
         truecp = twh_get_classical_potential()
         if ( truecp .eq. 'Repulsive Sphere' .and.
     &        localptype .eq. 'Hard Sphere' ) then
c              --- these use the same parameters
            lproblem = .false.
         endif
         if ( truecp .eq. 'Repulsive Well' .and.
     &        localptype .eq. 'Square Well' ) then
c              --- these use the same parameters
            lproblem = .false.
         endif
         if ( truecp .eq. 'Repulsive Multiwell' .and.
     &        localptype .eq. 'Multiwell' ) then
c              --- these use the same parameters
            lproblem = .false.
         endif
         if ( localptype .eq. 'wild' ) then
            lproblem = .false.
         endif
         if ( truecp .eq. 'Scaled Lennard-Jones' .and.
     &        localptype .eq. 'Lennard-Jones' ) then
c              --- these use the same parameters
            lproblem = .false.
         endif
         if ( truecp .eq. 'Weeks-Chandler-Anderson' 
     &        .and. localptype .eq. 'Lennard-Jones' ) then
c              --- these use the same parameters
            lproblem = .false.
         endif
         if ( truecp .eq. 'Repulsive 2580 Multistep' 
     &        .and. localptype .eq. 'Hard 2580 Multistep' ) then
            lproblem = .false.
         endif
         if ( truecp .eq. 'EAM pair only'
     &        .and. localptype .eq. 'Embedded Atom Method' ) then
            lproblem = .false.
         endif
         if ( truecp .eq. 'SW pair only'
     &        .and. localptype .eq. 'Stillinger-Weber' ) then
            lproblem = .false.
         endif
         if ( lproblem ) then
            write(6,*) 'RWFORCEFIELD: invalid potential mix'
            write(6,*) 'towhee_input: ',twh_get_classical_potential()
            write(6,*) 'ff_filename: ',localptype
            stop
         endif
      endif
   
      return
      end

      function twh_get_ncoeff(localptype, iversion, lrwhbond
     &     ,lgetcoeffnum, hbcoeffstart, hbcoeffstop)
c     ******************************************************************
c     * Return the number of coefficients for a given potential.       *
c     * In certain cases, lrwhbond, lgetcoeffnum, hbcoeffstart and     *
c     * hbcoeffstop are set as well.                                   *
c     *                                                                *
c     * split from rwforcefield 02-21-2006 by MAW                      *
c     * last modified 07-11-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NCOEFF
#include "functions.h"
c     --- variables passed to/from this subroutine
      integer iversion, hbcoeffstart, hbcoeffstop
      character*30 localptype
      logical lrwhbond, lgetcoeffnum

c     --- initialize
      twh_get_ncoeff = 0

      if ( localptype .eq. 'Lennard-Jones' ) then
c        --- Lennard-Jones 12-6 potential
c        --- first four nbcoeffs
         twh_get_ncoeff = 4
      elseif ( localptype .eq. '9-6' ) then
c        --- 9-6 Lennard Jones for Compass
c        --- 2 coeffs and also some extra information about
c        --- the atomnames
         twh_get_ncoeff = 2
      elseif ( localptype .eq. '12-6 plus 12-10 H-bond' ) then
c        --- essentially a LJ 12-6, but with some extra parameters in 
c        --- a separate data structure
         twh_get_ncoeff = 4
c        --- need the h-bond data structures
         lrwhbond = .true.
c        --- set the coeffcient boundaries for h-bond
         hbcoeffstart = 1
         hbcoeffstop = 3
      elseif ( localptype .eq. '12-6 plus solvation' ) then
c        --- no longer compatable with earlier versions with this 
c        --- localptype
         if ( iversion .le. 4 ) then
            write(6,*) 'RWFORCEFIELD: 12-6 + solvation not compatable'
            write(6,*) 'with earlier versions any longer'
            write(6,*) 'please update to version 5 of towhee_ff'
            stop
         endif
c        --- lj potential with an additional pairwise contribution 
c        --- in order to model water using implicit solvent 
c        --- model
         twh_get_ncoeff = 11
      elseif ( localptype .eq. '12-9-6' ) then
c        --- 12-9-6 potential
c        --- first three nbcoeffs
         twh_get_ncoeff = 3
      elseif ( localptype .eq. 'Buffered 14-7' ) then
c        --- first 4 nbcoeffs
         twh_get_ncoeff = 4
      elseif ( localptype .eq. 'Drukker H-bond' ) then
c        --- first 4 nbcoeffs
         twh_get_ncoeff = 4
      elseif ( localptype .eq. 'Double Exponential' ) then
c        --- first 4 nbcoeffs
         twh_get_ncoeff = 4
      elseif ( localptype .eq. 'Hard 2580 Multistep' .or.
     &        localptype .eq. 'Repulsive 2580 Multistep' ) then
c        --- 2580 potential, need the first 5 nbcoeffs
         twh_get_ncoeff = 6
      elseif ( localptype .eq. 'Embedded Atom Method' .or.
     &        localptype .eq. 'EAM pair only' ) then
c        --- embedded atom, these values are set later
         twh_get_ncoeff = 0
         lgetcoeffnum = .true.
      elseif ( localptype .eq. 'Exponential-6' ) then
c        --- exponential-6 potential 
c        --- must explicitly list all of the cross terms
         twh_get_ncoeff = 3
      elseif ( localptype .eq. 'Exponential-12-6' ) then
c        --- combination of exp-6 and repulsive 12
c        --- must explcitly declare cross terms
         twh_get_ncoeff = 4
      elseif ( localptype .eq. 'Gordon n-6' ) then
c        --- Peter Gordon's modified n-6 potential
         twh_get_ncoeff = 3
      elseif ( localptype .eq. 'Hard Sphere' .or. 
     &        localptype .eq. 'Repulsive Sphere' ) then
c        --- hard sphere
         twh_get_ncoeff = 1
      elseif ( localptype .eq. 'Multiwell' .or.
     &        localptype .eq. 'Repulsive Multiwell' ) then
c        --- stored in the tabulated data structures
         lgetcoeffnum = .true.
      elseif ( localptype .eq. 'Square Well' .or.
     &        localptype .eq. 'Repulsive Well' ) then
c        --- Square Well or repulsive well
         twh_get_ncoeff = 3
      elseif ( localptype .eq. 'Stillinger-Weber' .or.
     &        localptype .eq. 'SW pair only' ) then
c        --- stillinger-weber potential
         twh_get_ncoeff = 9
      elseif ( localptype .eq. 'Tabulated Pair' ) then
c        --- tabulated pair potential
         twh_get_ncoeff = 0
         lgetcoeffnum = .true.
      elseif ( localptype .eq. 'UFF 12-6' ) then
c        --- UFF Lennard-Jones potential with a bunch of other params
         twh_get_ncoeff = 12
      elseif ( localptype .eq. 'Weeks-Chandler-Anderson' ) then
c        --- truncated and shifted Lennard-Jones
         twh_get_ncoeff = 4
      elseif ( localptype .eq. 'wild' ) then
c        --- matches with everything, no coefficients
         twh_get_ncoeff = 0
      else
         write(6,*) 'RWFORCEFIELD: invalid classical_potential '
     &        ,localptype
         stop
      endif

      return
      end


      subroutine twh_get_vibration_params(iversion,type, vibstyle, 
     &   coeffstart, coeffstop)
c     ******************************************************************
c     * Returns the index limits for vibstyle coefficients, given a    *
c     * certain vibstyle.  type is used only for error output.         *
c     * Arguments modified by this subroutine:                         *
c     *  coeffstart                                                    *
c     *  coeffstop                                                     *
c     * Both arguments are initialized here.                           *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 07-18-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer iversion,type, vibstyle, coeffstart, coeffstop

      if ( vibstyle .eq. BOND_STYLE_FIXED ) then
c           --- fixed bond style
         coeffstart = 0
         coeffstop = 0
      elseif ( vibstyle .eq. BOND_STYLE_HARMONIC ) then
c           --- standard harmonic
         coeffstart = 0
         coeffstop = 1
      elseif ( vibstyle .eq. BOND_STYLE_GROMOS ) then
c           --- Gromos quartic
         coeffstart = 0
         coeffstop = 1
      elseif ( vibstyle .eq. BOND_STYLE_NONLINEAR ) then
c           --- nonlinear
         coeffstart = 0
         coeffstop = 2
      elseif ( vibstyle .eq. BOND_STYLE_MMTWO ) then
c           --- MM2
         coeffstart = 0
         coeffstop = 1
      elseif ( vibstyle .eq. BOND_STYLE_COMPASS ) then
c           --- Compass quartic
         coeffstart = 0
         coeffstop = 3
      elseif ( vibstyle .eq. BOND_STYLE_NONBOND ) then
c           --- pure nonbond
         coeffstart = 0
         coeffstop = -1
      elseif ( vibstyle .eq. BOND_STYLE_NONE ) then
c           --- no interactions
         coeffstart = 0
         coeffstop = -1
      elseif ( vibstyle .eq. BOND_STYLE_MORSE ) then
c           --- morse
         coeffstart = 0
         coeffstop = 2
      elseif ( vibstyle .eq. BOND_STYLE_SQUARE_WELL ) then
c        --- infinite square well
         if ( iversion .le. 14 ) then
c           --- had a dummy value of vib(0)
            coeffstart = 0
         else
c           --- vib(0) no longer required
            coeffstart = 1
         endif
         coeffstop = 3
      elseif ( vibstyle .eq. BOND_STYLE_HARMONIC_NONBOND ) then
c           --- standard harmonic plus nonbond
         coeffstart = 0
         coeffstop = 1
      elseif ( vibstyle .eq. BOND_STYLE_FENE ) then
c        --- FENE for CG models
         coeffstart = 1
         coeffstop = 2
      else
         write(6,*) 'RWFORCEFIELD: invalid vibstyle'
     &        ,type,vibstyle
         stop
      endif

      return
      end


      subroutine get_angle_params(anglestyle, iversion,
     &   coeffstart, coeffstop, ibastart, ibaend, ibbstart, ibbend,
     &   lrwbondangle, lbastring, lrwbondbond)
c     ******************************************************************
c     * obtain various parameters for given angle style.               *
c     * Arguments coeffstart, coeffstop are set always.                *
c     * Arguments ibastart, ibaend, ibbstart, ibbend,                  *
c     *  lrwbondangle, lbastring, lrwbondbond                          *
c     *  are set depending on angle style.                             *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer anglestyle, iversion, coeffstart, coeffstop
      integer ibastart, ibaend, ibbstart, ibbend
      logical lrwbondangle, lbastring, lrwbondbond

      lbastring = .false.

c     --- number of coefficients does not include the equilibrium
c     --- angle (coeff 0)
      if ( anglestyle .eq. ANG_STYLE_RIGID ) then
c           --- rigid angle
         coeffstart = 0
         coeffstop = 1
      elseif ( anglestyle .eq. ANG_STYLE_HARMONIC ) then
c           --- standard harmonic
         coeffstart = 0
         coeffstop = 1
      elseif ( anglestyle .eq. ANG_STYLE_DREIDING ) then
c           --- 1 + Cosine
         coeffstart = 0
         coeffstop = 1
      elseif ( anglestyle .eq. ANG_STYLE_HARM_COS ) then
c           --- harmonic cosine
         coeffstart = 0
         coeffstop = 1
      elseif ( anglestyle .eq. ANG_STYLE_COMPASS_AUTO ) then
c           --- Compass quartic with clever atom determination
         coeffstart = 0
         coeffstop = 3
c        --- output bond-angle and bond-bond information
         lrwbondangle = .true.
         ibastart = 4
         ibaend = 5
c        --- no longer used, but required to read the old versions
         lbastring = .true.
         lrwbondbond = .true.
         ibbstart = 6
         ibbend = 6
      elseif ( anglestyle .eq. ANG_STYLE_CHARMM_UB ) then
c           --- Charmm harmonic with Urey-Bradley
         coeffstart = 0
         coeffstop = 3
      elseif ( anglestyle .eq. ANG_STYLE_NONBOND ) then
c           --- nonbonded interactions, no coefficients
         coeffstart = 0
         coeffstop = -1
      elseif ( anglestyle .eq. ANG_STYLE_NB_HARM ) then
c           --- harmonic plus vdw and coulombic
         coeffstart = 0
         coeffstop = 1
      elseif ( anglestyle .eq. ANG_STYLE_COMPASS_EX ) then
c           --- Compass quartic with explicit terms
         coeffstart = 0
         coeffstop = 3
         lrwbondangle = .true.
         ibastart = 4
         ibaend = 7
         lrwbondbond = .true.
         ibbstart = 8
         ibbend = 10
      elseif ( anglestyle .eq. ANG_STYLE_FOURIER ) then
c           --- UFF style Fourier expansion plus single term
         coeffstart = 0
         coeffstop = 2
c           --- safety check as this style changed going from v9 to v10
         if ( iversion .lt. 10 ) then
            write(6,*) 'RWFORCEFIELD: anglestyle 9 was different'
     &           ,' in this old version and is not valid now'
            write(6,*) 'upgrade to a more current version of this'
     &           ,' forcefield'
            stop
         endif
      elseif ( anglestyle .eq. ANG_STYLE_TTFOURIER ) then
c           --- UFF style 3 term Fourier expansion
         coeffstart = 0
         coeffstop = 3
      elseif ( anglestyle .eq. ANG_STYLE_NONE ) then
c           --- no interaction
         coeffstart = 0
         coeffstop = -1
      elseif ( anglestyle .eq. ANG_STYLE_SIXTH_AUTO ) then
c           --- MM2 with cross terms
         coeffstart = 0
         coeffstop = 2
         ibastart = 3
         ibaend = 3
         lrwbondangle = .true.
c        --- no longer used, but required to read the old versions
         lbastring = .true.
      elseif ( anglestyle .eq. ANG_STYLE_SQUARE_WELL ) then
c           --- infinite square well
         coeffstart = 0
         coeffstop = 3
      elseif ( anglestyle .eq. ANG_STYLE_MULTI_RIGID ) then
c        --- multiple rigid angles
         coeffstart = 0
         coeffstop = 2
      elseif ( anglestyle .eq. ANG_STYLE_MMFF ) then
c        --- MMF cubic plus bond-stretch
         coeffstart = 0
         coeffstop = 2
         lrwbondangle = .true.
         ibastart = 3
         ibaend = 4
      elseif ( anglestyle .eq. ANG_STYLE_HARM_COS_OT ) then
c        --- harmonic cosine plus 1-3 nonbond
         coeffstart = 0
         coeffstop = 1
      else
         write(6,*) 'RWFORCEFIELD: invalid anglestyle'
     &        ,anglestyle
         stop
      endif
      return
      end


      subroutine twh_get_torsion_params(lfailure,iversion, rwflag, 
     &   type, coeffstart, coeffstop, ltorstring)
c     ******************************************************************
c     * obtain various parameters for given torsion style.             *
c     * Arguments set in this procedure:                               *
c     *  coeffstart                                                    *
c     *  coeffend                                                      *
c     *  ltorstring                                                    *
c     * (values of these arguments initialized here)                   *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NTORLOOP
#define FUNCTION_GET_TORSTYLE
#include "functions.h"
c     --- variables passed to/from this subroutine
      logical lfailure,ltorstring
      integer iversion, rwflag, type
      integer coeffstart, coeffstop

c     --- local variables
      character*15 idname
      integer torstyle

      idname = 'rwforcefield'
      ltorstring = .false.
      torstyle = twh_get_torstyle(type)
      if ( torstyle .eq. 1 ) then
c           --- simple harmonic
         coeffstart = 0
         if ( iversion .lt. 11 ) then
c              --- there was a mistake in this old version that had
c              --- an extra coeff that was not used
            coeffstop = 2
         else
            coeffstop = 1
         endif
      elseif ( torstyle .eq. 2 ) then
c           --- OPLS style cosine series
         coeffstart = 1
         coeffstop = 3
      elseif ( torstyle .eq. 3 ) then
c           --- Gromos/Charmm/Amber cosine series
         call twh_rwtorloop(lfailure,type,rwflag)
         if ( lfailure ) return
         coeffstart = 1
         coeffstop = 3*twh_get_ntorloop(type)
      elseif ( torstyle .eq. 4 ) then
c           --- Gromos/Charmm/Amber cosine series
c           --- with additional harmonic
         call twh_rwtorloop(lfailure,type,rwflag)
         if ( lfailure ) return
         coeffstart = 1
         coeffstop = 3*twh_get_ntorloop(type)+2
      elseif ( torstyle .eq. 5 ) then
c           --- compass with all of the cross terms
         coeffstart = 0
         coeffstop = 19
c           --- torstrings
         ltorstring = .true.
      elseif ( torstyle .eq. 6 ) then
c           --- compass without the cross terms
         coeffstart = 0
         coeffstop = 2
      elseif ( torstyle .eq. 7 ) then
c           --- TraPPE simple cosine function
         coeffstart = 0
         coeffstop = 1
      elseif ( torstyle .eq. 8 ) then
c           --- nonbonded potential only
         coeffstart = 0
         if ( iversion .lt. 12 ) then
c           --- there was an unused coefficient in previous version
            coeffstop = 0
         else
c           --- do not need anything in addition to the 1-4 scaling
            coeffstop = -1
         endif
      elseif ( torstyle .eq. 9 ) then
c           --- class 2 with explicit terms
         coeffstart = 0
         coeffstop = 31
      elseif ( torstyle .eq. 10 ) then
c           --- cosine power series
         call twh_rwtorloop(lfailure,type,rwflag)
         if ( lfailure ) return
         coeffstart = 0
         coeffstop = twh_get_ntorloop(type)
      elseif ( torstyle .eq. 11 ) then
c           --- old UA OPLS style cosine series
         coeffstart = 0
         coeffstop = 3
      elseif ( torstyle .eq. 12 ) then
c           --- NERD Version 3 style cosine series
         call twh_rwtorloop(lfailure,type,rwflag)
         if ( lfailure ) return
         coeffstart = 0
         coeffstop = twh_get_ntorloop(type)
      elseif ( torstyle .eq. 13 ) then
c           --- two term cosine series from OPLS-ua amides
         coeffstart = 1
         coeffstop = 2
      elseif ( torstyle .eq. 14 ) then
c           --- UFF cosine
         coeffstart = 1
         coeffstop = 3
      elseif ( torstyle .eq. 15 ) then
c           --- DREIDING cosine
         coeffstart = 1
         coeffstop = 3
      elseif ( torstyle .eq. 16 ) then
c           --- 2-fold cosine
         coeffstart = 1
         coeffstop = 1
      elseif ( torstyle .eq. 17 ) then
c           --- potoff cosine
         coeffstart = 1
         coeffstop = 3
      elseif ( torstyle .eq. 18 ) then
c           --- square well
         coeffstart = 1
         coeffstop = 3
      elseif ( torstyle .eq. 19 ) then
c           --- Amber cosine series divided by number of torsions
         call twh_rwtorloop(lfailure,type,rwflag)
         if ( lfailure ) return
         coeffstart = 1
         coeffstop = 3*twh_get_ntorloop(type)
      elseif ( torstyle .eq. 20 ) then
c        --- OPLS 4 parameter style cosine series
         coeffstart = 1
         coeffstop = 4
      elseif ( torstyle .eq. 21 ) then
         call twh_rwtorloop(lfailure,type,rwflag)
         if ( lfailure ) return
         coeffstart = 0
         coeffstop = twh_get_ntorloop(type)
      elseif ( torstyle .eq. 22 ) then
         coeffstart = 0
         coeffstop = 5
      elseif ( torstyle .eq. 23 ) then
         coeffstart = 0
         coeffstop = 1
      else
         write(6,*) 'RWFORCEFIELD: invalid torstyle'
     &        ,type,torstyle
         stop
      endif
 
      return
      end


      subroutine get_improper_params(impstyle, coeffstart, coeffstop)
c     ******************************************************************
c     * obtain various parameters for given improper torsion style.    *
c     * Arguments set in this procedure:                               *
c     *  coeffstart                                                    *
c     *  coeffstop                                                     *
c     * (values of these arguments initialized here)                   *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 12-17-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer impstyle, coeffstart, coeffstop

      if ( impstyle .eq. 1 ) then
c           --- simple harmonic
         coeffstart = 0
         coeffstop = 1
      elseif ( impstyle .eq. 2 ) then
c           --- OPLS style cosine series
         coeffstart = 1
         coeffstop = 3
      elseif ( impstyle .eq. 3 ) then
c           --- MGM sterochem inforcer
         coeffstart = 0
         coeffstop = 1
      elseif ( impstyle .eq. 4 ) then
c           --- Amber improper
         coeffstart = 1
         coeffstop = 3
      elseif ( impstyle .eq. 5 ) then
c           --- Wilson harmonic
         coeffstart = 0
         coeffstop = 1
      elseif ( impstyle .eq. 6 ) then
c           --- UFF cosine series
         coeffstart = 0
         coeffstop = 3
      elseif ( impstyle .eq. 7 ) then
c           --- 1-cosine
         coeffstart = 1
         coeffstop = 1
      elseif ( impstyle .eq. 8 ) then
         coeffstart = 1
         coeffstop = 1
      else
         write(6,*) 'RWFORCEFIELD: invalid impstyle' ,impstyle
         stop
      endif
   
      return
      end

      subroutine twh_rw_eam_params(lfailure,rwflag,type,ntypes
     &     ,iversion)
c     ******************************************************************
c     * obtain various parameters for given improper torsion style.    *
c     * Arguments set in this procedure:                               *
c     * Global arguments set in this procedure:                        *
c     *                                                                *
c     * split from rwforcefield 02-22-2006 by MAW                      *
c     * last modified 07-16-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag, type, ntypes, iversion
c     --- local variables
      logical lneedthree
      character*15 idname
      character*20 eamembedstyle,eamdensstyle
      character*50 goodstring
      integer jtype,itest,jtest,iscan,eam_ndens,eam_nembed
      double precision avalue,bvalue,cvalue

      idname = 'rwforcefield'
      if ( iversion .lt. 9 ) then
c        --- old format for the embedding density depended only upon
c        --- the contributing type.  Only need to read this as we no
c        --- longer write older versions
         if ( rwflag .eq. WRITE_FLAG ) then
            write(6,*) 'RWFORCEFIELD: writing old version'
            stop
         endif
         goodstring = 'eam_dens'
         call twh_read_label(lfailure,FFUNIT,goodstring, .false.,idname)
         if ( lfailure ) return
         read(FFUNIT,*) itest,eam_ndens
         call twh_eam_ndens(GLB_SET,type,type,eam_ndens)
c        --- safety check
         if ( itest .ne. type ) then
            write(6,*) 'RWEAM: mismatch of eam_dens types'
            write(6,*) 'type,itest',type,itest
            stop
         endif
         if ( eam_ndens .gt. MAXTABLE ) then
            write(6,*) 'RWEAM: eam_ndens > MAXTABLE'
            write(6,*) 'type,eam_ndens,MAXTABLE'
     &            ,type,eam_ndens,MAXTABLE
            stop
         endif

c        --- eam_dens_style
         goodstring = 'eam_dens_style'
         call twh_read_labeled_string(lfailure,FFUNIT, goodstring
     &        , .false.,idname,eamdensstyle)
         if ( lfailure ) return
         call twh_eam_dens_style(GLB_SET,type,type,eamdensstyle)

c        --- eam_dens data
         goodstring = 'eam_dens_data'
         call twh_read_label(lfailure,FFUNIT,goodstring, .false.,idname)
         if ( lfailure ) return
         call twh_eam_ndens(GLB_GET,type,type,eam_ndens)
         do iscan = 1,eam_ndens
            read(FFUNIT,*) avalue,bvalue
            call twh_eam_dens(GLB_SET,1,type,type,iscan,avalue)
            call twh_eam_dens(GLB_SET,2,type,type,iscan,bvalue)
         enddo
c        --- set all of the cross densities to this density
         call setcrossdens(type)
      else
c        --- version 9 and later has a full cross term table for the
c        --- embedding density
         do jtype = 1,ntypes
            goodstring = 'eam_ndens'
            if ( rwflag .eq. WRITE_FLAG ) then
               write(FFUNIT,'(a9)') goodstring
               call twh_eam_ndens(GLB_GET,jtype,type,eam_ndens)
               write(FFUNIT,*) jtype,type,eam_ndens
            else
               call twh_read_label(lfailure,FFUNIT,goodstring
     &              , .false.,idname)
               if ( lfailure ) return
               read(FFUNIT,*) jtest,itest,eam_ndens
               call twh_eam_ndens(GLB_SET,jtype,type,eam_ndens)
c              --- safety check
               if ( itest .ne. type .or. jtest .ne. jtype ) then
                  write(6,*) 'RWFORCEFIELD: mismatch of density'
     &                  ,' types'
                  write(6,*) 'jtest,jtype',jtest,jtype
                  write(6,*) 'itest,type',itest,type
                  stop
               endif
            endif
            if ( eam_ndens .gt. MAXTABLE ) then
               write(6,*) 'RWFORCEFIELD: eam_ndens > MAXTABLE'
               write(6,*) 'jtype,type,eam_ndens,MAXTABLE'
     &               ,jtype,type,eam_ndens,MAXTABLE
               stop
            endif

c           --- eam_dens_style
            goodstring = 'eam_dens_style'
            if ( rwflag .eq. WRITE_FLAG ) then
               write(FFUNIT,'(a14)') goodstring
               call twh_eam_dens_style(GLB_GET,jtype,type,eamdensstyle)
               write(FFUNIT,'(a20)') eamdensstyle
            else
               call twh_read_labeled_string(lfailure, FFUNIT
     &              , goodstring,.false.,idname,eamdensstyle)
               if ( lfailure ) return
               call twh_eam_dens_style(GLB_SET,jtype,type,eamdensstyle)
            endif

c           --- eam_dens data
            goodstring = 'eam_dens_data'
            if ( rwflag .eq. WRITE_FLAG ) then
               write(FFUNIT,'(a13)') goodstring
            else
               call twh_read_label(lfailure,FFUNIT,goodstring
     &              , .false.,idname)
               if ( lfailure ) return
            endif
            call twh_eam_ndens(GLB_GET,jtype,type,eam_ndens)
            do iscan = 1,eam_ndens
               if ( rwflag .eq. WRITE_FLAG ) then
                  call twh_eam_dens(GLB_GET,1,jtype,type,iscan,avalue)
                  call twh_eam_dens(GLB_GET,2,jtype,type,iscan,bvalue)
                  write(FFUNIT,'(e20.10,1x,e20.10)') avalue,bvalue
               else
                  read(FFUNIT,*) avalue,bvalue
                  call twh_eam_dens(GLB_SET,1,jtype,type,iscan,avalue)
                  call twh_eam_dens(GLB_SET,2,jtype,type,iscan,bvalue)
               endif
            enddo
         enddo
      endif

c     --- eam_embed nembed
      goodstring = 'eam_embed'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a9)') goodstring
         call twh_eam_nembed(GLB_GET,type,eam_nembed)
         write(FFUNIT,*) type,eam_nembed
      else
         call twh_read_label(lfailure,FFUNIT,goodstring, .false.,idname)
         if ( lfailure ) return
         read(FFUNIT,*) itest,eam_nembed
         call twh_eam_nembed(GLB_SET,type,eam_nembed)
c              --- safety check
         if ( itest .ne. type ) then
            write(6,*) 'RWEAM: mismatch of eam_embed types'
            write(6,*) 'type,itest',type,itest
            stop
         endif
         if ( eam_nembed .gt. MAXTABLE ) then
            write(6,*) 'RWEAM: eam_nembed > MAXTABLE'
            write(6,*) 'type,eam_nembed,MAXTABLE'
     &            ,type,eam_nembed,MAXTABLE
            stop
         endif
      endif

c     --- eam_embed_style
      goodstring = 'eam_embed_style'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a15)') goodstring
         call twh_eam_embed_style(GLB_GET,type,eamembedstyle)
         write(FFUNIT,'(a20)') eamembedstyle
      else
         call twh_read_labeled_string(lfailure, FFUNIT,
     &        goodstring, .false.,idname,eamembedstyle)
         if ( lfailure ) return
         call twh_eam_embed_style(GLB_SET,type,eamembedstyle)
      endif

c     --- check eam_embed_style to see if it one that requires 3
      if ( eamembedstyle .eq. 'Belashchenko 10 rho' ) then
         lneedthree = .true.
      else
         lneedthree = .false.
      endif

c     --- eam_embed data
      goodstring = 'eam_embed_data'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(FFUNIT,'(a14)') goodstring
      else
         call twh_read_label(lfailure,FFUNIT,goodstring, .false.,idname)
         if ( lfailure ) return
      endif
      call twh_eam_nembed(GLB_GET,type,eam_nembed)
      do iscan = 1,eam_nembed
         if ( rwflag .eq. WRITE_FLAG ) then
            call twh_eam_embed(GLB_GET,1,type,iscan,avalue)
            call twh_eam_embed(GLB_GET,2,type,iscan,bvalue)
            if ( lneedthree ) then
               call twh_eam_embed(GLB_GET,3,type,iscan,cvalue)
               write(FFUNIT,'(e20.10,1x,e20.10,1x,e20.10)')
     &              avalue,bvalue,cvalue
            else
               write(FFUNIT,'(e20.10,1x,e20.10)') avalue,bvalue
            endif
         else
            if ( lneedthree ) then
               read(FFUNIT,*) avalue,bvalue,cvalue
            else
               read(FFUNIT,*) avalue,bvalue
            endif
            call twh_eam_embed(GLB_SET,1,type,iscan,avalue)
            call twh_eam_embed(GLB_SET,2,type,iscan,bvalue)
            if ( lneedthree ) then
               call twh_eam_embed(GLB_SET,3,type,iscan,cvalue)
            endif
         endif
      enddo

      return
      end
