#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2014 Marcus G. Martin                       *
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
      subroutine twh_writedlpoly(lfailure)
c     ******************************************************************
c     * creates input files CONFIG and FIELD for DL_POLY, a molecular  *
c     * dynamics simulation code available from Daresbury labs, UK.    *
c     * The Dl_POLY code can be obtained at                            *
c     * http://www.cse.clrc.ac.uk/msi/software/DL_POLY                 *
c     * Note that user still needs to prepare the CONTROL file         *
c     *                                                                *
c     * originally written 08-05-2006 by O. Yazaydin                   *
c     * modified 12-17-2007 by C.M. Tenney                             *
c     *     angle types 2 (DREIDING) and 3 (harmonic cosine)           *
c     *     dihedral types 15 (DREIDING) and (not) 8 (nonbond)         *
c     *     improper form 4 type 7 (DREIDING 1-cos)                    *
c     *     fixed unit cell vector output and changed to triclinic pbc *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_ANGLESTYLE
#define FUNCTION_GET_BENCOEFF
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_ITBEN
#define FUNCTION_GET_ITTOR
#define FUNCTION_GET_LHERE
#define FUNCTION_GET_MASS
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_GET_NBOXI
#define FUNCTION_GET_NCMT
#define FUNCTION_GET_TORCOEFF
#define FUNCTION_GET_TORSTYLE
#define FUNCTION_GET_VIBCOEFF
#define FUNCTION_OLDGETATOMNUM
#define FUNCTION_GET_MASS
#define FUNCTION_GET_NBNAME
#define FUNCTION_GET_IMPCOEFF
#define FUNCTION_GET_IMPSTYLE
#define FUNCTION_GET_ITIMPROP
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
      logical ladd,lconstraints
      integer dlp_bond,dlp_angle,dlp_dihedral
      integer nuni,testmolty,reset,itype,jtype,ij
      integer ibox,ivib,ibend,itor,type,ichain,imolty,iunit
      integer count_atom,count_bond,nchain,nunit
      integer style,inpnum,count_angle,count_dihedral,nmolty
      integer dlp_improper, count_improper, ncmt, icmt,numboxes
      character*30 filename,localcp
      double precision xtemp,ytemp,ztemp,radtodeg,charge
      double precision totbond,hm1,hm2,hm3

c     --- retrieve constants
      call twh_constant_radtodeg(GLB_GET,radtodeg)
      call twh_nchain(GLB_GET,nchain)
      call twh_numboxes(GLB_GET,numboxes)

      do ibox=1, numboxes
c      --- start of FIELD file
c      --- open FIELD
       inpnum = 23
       write (filename, '(a6,i1)') 'FIELD_',ibox
       call twh_openfile(lfailure,inpnum,filename,'DEFAULT','UNKNOWN'
     &      ,'DEFAULT')
       if ( lfailure ) return
       write(inpnum,'(a36)') 'DL_POLY FIELD file written by Towhee'
       write(inpnum,'(a8)') 'UNITS kj'
       call twh_nmolty(GLB_GET,nmolty)
       write(inpnum,'(a10,i2)') 'MOLECULES ', nmolty
       do imolty=1,nmolty
          call twh_nunit(GLB_GET,imolty,nunit)
          write(inpnum,'(a17,i2)') 'Name of Molecule ', imolty
          write(inpnum,'(a9,i6)') 'NUMMOLS  '
     &         , twh_get_ncmt(ibox,imolty)
          write(inpnum,'(a7,i6)') 'ATOMS  ', nunit
          do iunit = 1,nunit
             call twh_ntype(GLB_GET,imolty,iunit,type)
             call twh_qqatom(GLB_GET,imolty,iunit,charge)
             write(inpnum,'(a8,2f12.5,2a5)')
     &            twh_get_nbname(type),twh_get_mass(type)
     &            ,charge,'    1','    0'
          enddo
c         ---- number of bonds,angles,torsions
          dlp_bond = 0
          dlp_angle = 0
          dlp_dihedral = 0
          dlp_improper = 0
          do ichain = 1,nchain
             if ( twh_get_nboxi(ichain) .eq. ibox ) then
                call twh_moltyp(GLB_GET,ichain,testmolty)
                if ( imolty .eq. testmolty ) then
                   do iunit = 1,nunit
c                     --- number of bonds
                      do ivib = 1,twh_get_invib(imolty,iunit)
                         if ( twh_get_ijvib(imolty,iunit,ivib) .gt.
     &                        iunit ) then
                            dlp_bond = dlp_bond + 1
                         endif
                      enddo
c                     --- number of angles
                      do ibend = 1,twh_get_inben(imolty,iunit)
                         type = twh_get_itben(imolty,iunit,ibend)
                         ladd = .false.
c                        --- add based solely on ijben3 > iunit
                         if ( twh_get_ijben3(imolty,iunit,ibend)
     &                        .gt. iunit ) ladd = .true.
                         if ( ladd ) then
                            dlp_angle = dlp_angle + 1
                         endif
                      enddo
c                     --- number of torsions
                      do itor = 1,twh_get_intor(imolty,iunit)
                         type = twh_get_ittor(imolty,iunit,itor)
                         if ( type .gt. 0 ) then
                            ladd = .false.
c                           --- any other case compare ijt4 > iunit
                            if ( twh_get_ijtor4(imolty,iunit,itor)
     &                           .gt. iunit ) then
                               style = twh_get_torstyle(type)
                               if ( style .ne. 8) then
                                  ladd = .true.
                               endif
                            endif
                            if ( ladd ) then
                               dlp_dihedral = dlp_dihedral + 1
                            endif
                         endif
                      enddo
c                     --- end number of torsions

c                     --- number of impropers
                      do itor = 1,twh_get_inimprop(imolty,iunit)
                         type = twh_get_itimprop(imolty,iunit,itor)
                         if ( type .gt. 0 ) then
c c                         --- compare ijimprop4 > iunit
                            dlp_improper = dlp_improper + 1
                         endif
                      enddo
c                     --- end number of impropers
                   enddo
                   goto 111
                endif
             endif
          enddo
 111      continue
c         ---  output bond information
          if ( dlp_bond .ne. 0 ) then
             count_bond = 0
c            --- search for CONSTRAINTS
c            --- DL_POLY treats fixed bonds with SHAKE
c            --- and lists under a seperate section called CONSTRAINTS
             lconstraints = .false.
             do ichain = 1,nchain
                if ( twh_get_nboxi(ichain) .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      reset = twh_oldgetatomnum(ichain,1)-1
                      do iunit = 1,nunit
                         do ivib = 1,twh_get_invib(imolty,iunit)
                            call twh_itvib(GLB_GET,imolty,iunit,ivib
     &                           ,type)
                            if (type .gt.0) then
                               if ( twh_get_ijvib(imolty,iunit,ivib)
     &                              .gt. iunit ) then
                                  call twh_vibstyle(GLB_GET,type
     &                                 ,style)
                                  if ( style .eq. BOND_STYLE_FIXED) then
                                     count_bond = count_bond + 1
                                     if (count_bond .eq. 1) then
                                        lconstraints = .true.
                                        write(inpnum,'(a13,i5)')
     &                                       'CONSTRAINTS  ', dlp_bond
                                     endif
                                     write(inpnum,'(2i5,f12.5)')
     &                                    twh_oldgetatomnum(ichain
     &                                    ,iunit) - reset
     &                                    ,twh_oldgetatomnum(ichain
     &                                    ,twh_get_ijvib(imolty,iunit
     &                                    ,ivib)) - reset
     &                                    ,twh_get_vibcoeff(type,0)
                                  endif
                               endif
                            endif
                         enddo
                      enddo
                      goto 222
                   endif
                endif
             enddo
 222         continue
c            --- look for flexible bonds
c            --- subtract # of fixed bonds from flexible bonds
             if (count_bond .gt. 0) dlp_bond = dlp_bond - count_bond
c            --- if there are no flexible bonds continue with angles
             if (dlp_bond .eq. 0) goto 444
             write(inpnum,'(a8,i5)') 'BONDS   ', dlp_bond
             do ichain = 1,nchain
                if ( twh_get_nboxi(ichain) .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      reset = twh_oldgetatomnum(ichain,1)-1
                      do iunit = 1,nunit
                         do ivib = 1,twh_get_invib(imolty,iunit)
                            call twh_itvib(GLB_GET,imolty,iunit,ivib
     &                           ,type)
                            if (type .gt.0) then
                               if ( twh_get_ijvib(imolty,iunit,ivib)
     &                              .gt. iunit ) then
                                  count_bond = count_bond + 1
                                  call twh_vibstyle(GLB_GET,type
     &                                 ,style)
                                  if ( style .eq. BOND_STYLE_HARMONIC
     &                                 ) then
                                     write(inpnum,'(a5,2i5,2f12.5)')
     &                                    'harm '
     &                                    ,twh_oldgetatomnum(ichain
     &                                    ,iunit) - reset
     &                                    ,twh_oldgetatomnum(ichain
     &                                    ,twh_get_ijvib(imolty,iunit
     &                                    ,ivib)) - reset
     &                                    ,twh_get_vibcoeff(type,1)
     &                                    *2.0d0/CON_KJTOK
     &                                    ,twh_get_vibcoeff(type,0)
                                  else
                                     if ( .not. lconstraints ) then
                                        write(6,*) 'WRITEDLPOLY: ',
     &                                       'invalid bond style '
     &                                       ,style,' has nothing '
     &                                       ,'similar in dlpoly'
                                     endif
                                  endif
                               endif
                            endif
                         enddo
                      enddo
                      goto 333
                   endif
                endif
             enddo
 333         continue
             if ( dlp_bond .ne. count_bond ) then
                write(6,*) 'Error in writedlpoly'
                write(6,*) 'dlp_bond ',dlp_bond
                write(6,*) 'count_bond ',count_bond
             endif
          endif
 444      continue
c         ----  output angle information
          write(inpnum,'(a8,i5)') 'ANGLES   ', dlp_angle
          if ( dlp_angle .ne. 0 ) then
             count_angle = 0
             do ichain = 1,nchain
                if ( twh_get_nboxi(ichain) .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      reset = twh_oldgetatomnum(ichain,1)-1
                      do iunit = 1,nunit
                         do ibend = 1,twh_get_inben(imolty,iunit)
                            type = twh_get_itben(imolty,iunit,ibend)
                            if ( twh_get_ijben3(imolty,iunit,ibend)
     &                           .gt.iunit ) then
                               count_angle = count_angle + 1
                               style = twh_get_anglestyle(type)
                               if ( style .eq. ANG_STYLE_HARMONIC ) then
                                  write(inpnum,'(a5,3i5,2f12.5)')
     &                                 'harm '
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,iunit) - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,twh_get_ijben2(imolty,iunit
     &                                 ,ibend))
     &                                 - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,twh_get_ijben3(imolty,iunit
     &                                 ,ibend))
     &                                 - reset
     &                                 ,twh_get_bencoeff(type,1)
     &                                 *2.0d0/CON_KJTOK
     &                                 ,twh_get_bencoeff(type,0)
     &                                 *radtodeg
                               elseif ( style .eq. ANG_STYLE_DREIDING
     &                                 ) then
                                  write(inpnum,'(a5,3i5,3f12.5)')
     &                                 'cos  '
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,iunit) - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,twh_get_ijben2(imolty,iunit
     &                                 ,ibend))
     &                                 - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,twh_get_ijben3(imolty,iunit
     &                                 ,ibend))
     &                                 - reset
     &                                 ,twh_get_bencoeff(type,1)
     &                                 *2.0d0/CON_KJTOK
     &                                 , 0.0, 1.0
                               elseif ( style .eq. ANG_STYLE_HARM_COS
     &                                 ) then
                                  write(inpnum,'(a5,3i5,2f12.5)')
     &                                 'hcos '
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,iunit) - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,twh_get_ijben2(imolty,iunit
     &                                 ,ibend))
     &                                 - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,twh_get_ijben3(imolty,iunit
     &                                 ,ibend))
     &                                 - reset
     &                                 ,twh_get_bencoeff(type,1)
     &                                 *2.0d0/CON_KJTOK
     &                                 ,twh_get_bencoeff(type,0)
     &                                 *radtodeg
                               else
                                  write(6,*) 'WRITEDLPOLY: ',
     &                                 'invalid angle style ',style, 
     &                                 ' has nothing similar '
     &                                 ,'in dlpoly'
                               endif
                            endif
                         enddo
                      enddo
                      goto 555
                   endif
                endif
             end do
 555         continue
             if ( dlp_angle .ne. count_angle ) then
                write(6,*) 'problem with angles in writelammps'
                write(6,*) 'dlp_angle ',dlp_angle
                write(6,*) 'count_angle ',count_angle
             endif
          endif
c         ---- output torsion information
          write(inpnum,'(a11,i5)') 'DIHEDRALS  ', dlp_dihedral
          if ( dlp_dihedral .ne. 0 ) then
             count_dihedral = 0
             do ichain = 1,nchain
                if ( twh_get_nboxi(ichain) .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      reset = twh_oldgetatomnum(ichain,1)-1
                      do iunit = 1,nunit
                         do itor = 1,twh_get_intor(imolty,iunit)
                            type = twh_get_ittor(imolty,iunit,itor)
                            if ( type .gt. 0 ) then
                               if ( twh_get_ijtor4(imolty,iunit,itor)
     &                              .gt. iunit ) then
c                                 --- ignore style 8 as it has
c                                 --- (one-four interactions only)
                                  style = twh_get_torstyle(type)
                                  if ( style .ne. 8) then
                                     count_dihedral = count_dihedral
     &                                    + 1
                                     if ( style .eq. 2) then
                                        write(inpnum
     &                                       ,'(a5,4i5,3f12.5,2a8)')
     &                                       'cos3 '
     &                                       ,twh_oldgetatomnum(
     &                                       ichain,iunit) - reset
     &                                       ,twh_oldgetatomnum(
     &                                       ichain
     &                                       ,twh_get_ijtor2(imolty
     &                                       ,iunit,itor)) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,twh_get_ijtor3(imolty
     &                                       ,iunit,itor)) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,twh_get_ijtor4(imolty
     &                                       ,iunit,itor)) - reset
     &                                       ,twh_get_torcoeff(type,1)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,twh_get_torcoeff(type,2)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,twh_get_torcoeff(type,3)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,'     0.5', '     0.5'
                                     else if ( style .eq. 11) then
                                        write(inpnum
     &                                       ,'(a5,4i5,3f12.5,2a8
     &,f12.5)') 'opls '
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,iunit) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,twh_get_ijtor2(imolty
     &                                       ,iunit,itor))
     &                                       - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,twh_get_ijtor3(imolty
     &                                       ,iunit,itor)) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,twh_get_ijtor4(imolty
     &                                       ,iunit,itor)) - reset
     &                                       ,twh_get_torcoeff(type,0)
     &                                       /CON_KJTOK
     &                                       ,twh_get_torcoeff(type,1)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,twh_get_torcoeff(type,2)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,'     0.5', '     0.5'
     &                                       ,twh_get_torcoeff(type,3)
     &                                       *2.0d0/CON_KJTOK
                                     else if ( style .eq. 15) then
                                        totbond =
     &                                       dble((twh_get_invib(
     &                                       imolty,twh_get_ijtor2(
     &                                       imolty,iunit,itor))-1)
     &                                       *(twh_get_invib(imolty
     &                                       ,twh_get_ijtor3(imolty
     &                                       ,iunit,itor))-1))
                                        write(inpnum,'(a5,4i5,3f12.5,2a8
     &,f12.5)')
     &                                       'cos '
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,iunit) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,twh_get_ijtor2(imolty
     &                                       ,iunit,itor)) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,twh_get_ijtor3(imolty
     &                                       ,iunit,itor))
     &                                       - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,twh_get_ijtor4(imolty
     &                                       ,iunit,itor)) - reset
     &                                       ,twh_get_torcoeff(type,1)
     &                                       /CON_KJTOK/totbond
     &                                       ,twh_get_torcoeff(type,2)
     &                                       *twh_get_torcoeff(type,3)
     &                                       *radtodeg-180.0d0
     &                                       ,twh_get_torcoeff(type,2)
     &                                       ,'     1.0', '     1.0'
                                     else
                                        write(6,*) 'WRITEDLPOLY: ',
     &                                       'invalid torsion style '
     &                                       ,style,
     &                                       ' has nothing similar'
     &                                       ,' in dlpoly'
                                     endif
                                  endif
                               endif
                            endif
                         enddo
                      enddo
                      goto 666
                   endif
                endif
             enddo
 666         continue
             if ( dlp_dihedral .ne. count_dihedral ) then
                write(6,*) 'problem with angles in writedlpoly'
                write(6,*) 'dlp_dihedral ',dlp_dihedral
                write(6,*) 'count_dihedral ',count_dihedral
             endif
          endif

c         ---- output improper torsion (inversion) information
          write(inpnum,'(a11,i5)') 'INVERSIONS  ', dlp_improper
          if ( dlp_improper .ne. 0 ) then
             count_improper = 0
             do ichain = 1,nchain
                if ( twh_get_nboxi(ichain) .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      reset = twh_oldgetatomnum(ichain,1)-1
                      do iunit = 1,nunit
                         do itor = 1,twh_get_inimprop(imolty,iunit)
                            type = twh_get_itimprop(imolty,iunit,itor)
                            if ( type .gt. 0 ) then
                               count_improper = count_improper + 1
                               style = twh_get_impstyle(type)
                               if (style .eq. 7) then
                                  write(inpnum,'(a5,4i5,f12.5)')
     &                                 'plan '
     &                                 ,twh_oldgetatomnum(ichain,iunit)
     &                                 - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,twh_get_ijimprop2(imolty,iunit
     &                                 ,itor))
     &                                 - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,twh_get_ijimprop3(imolty,iunit
     &                                 ,itor))
     &                                 - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,twh_get_ijimprop4(imolty,iunit
     &                                 ,itor))
     &                                 - reset
     &                                 ,twh_get_impcoeff(type,1)
     &                                 /CON_KJTOK
                               else
                                  write(6,*) 'WRITEDLPOLY: ',
     &                                 'invalid improper style ',style,
     &                                 ' has nothing similar in dlpoly'
                               endif
                            endif
                         enddo
                      enddo
                      goto 777
                   endif
                endif
             enddo
 777         continue
             if ( dlp_improper .ne. count_improper ) then
                write(6,*) 'problem with angles in writedlpoly'
                write(6,*) 'dlp_improper ',dlp_improper
                write(6,*) 'count_improper ',count_improper
             endif
          endif

          write(inpnum,'(a6)') 'FINISH'
       enddo
c      --- ouptput nonbonded parameters
c      --- get the classical potential
       call twh_classical_potential(GLB_GET,localcp)
c      --- determine the number of unique nonbonded interactions (nuni)
       nuni=0
       do itype = 1,NNTYPE
          do jtype = itype,NNTYPE
             if ( twh_get_lhere(itype)
     &            .and. twh_get_lhere(jtype) ) then
                nuni=nuni+1
             endif
          enddo
       enddo
       write(inpnum,'(a6,i5)') 'VDW   ', nuni
       if ( localcp .eq. 'Lennard-Jones' ) then
          do itype = 1,NNTYPE
             do jtype = itype,NNTYPE
                if ( twh_get_lhere(itype)
     &               .and. twh_get_lhere(jtype) ) then
                   ij = twh_getnbtype(itype,jtype)
                   write(inpnum,'(2a8,1x,a4,2f12.5)')
     &                  twh_get_nbname(itype)
     &                  ,twh_get_nbname(jtype),'lj  '
     &                  ,twh_get_nbcoeff(14,ij)/(4.0d0*CON_KJTOK)
     &                  ,(twh_get_nbcoeff(13,ij))**(1.0/6.0)
                endif
             enddo
          enddo
       end if
       write(inpnum,'(a5)') 'CLOSE'
       close (inpnum)
      enddo
c     --- end of FIELD file
c     --- start of CONFIG file
      do ibox=1,numboxes
         inpnum = 33
         write (filename, '(a7,i1)') 'CONFIG_',ibox
c        --- filename = 'CONFIG'
         call twh_openfile(lfailure,inpnum,filename,'DEFAULT','UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) return
         write(inpnum,'(a37)') 'DL_POLY CONFIG file written by Towhee'
         write(inpnum,'(2a10)') '         0','         3'
c        --- output box dimensions
         call twh_hmatrix(GLB_GET,ibox,1,1,hm1)
         call twh_hmatrix(GLB_GET,ibox,1,2,hm2)
         call twh_hmatrix(GLB_GET,ibox,1,3,hm3)
         write(inpnum,'(3f20.12)') hm1,hm2,hm3
         call twh_hmatrix(GLB_GET,ibox,2,1,hm1)
         call twh_hmatrix(GLB_GET,ibox,2,2,hm2)
         call twh_hmatrix(GLB_GET,ibox,2,3,hm3)
         write(inpnum,'(3f20.12)') hm1,hm2,hm3
         call twh_hmatrix(GLB_GET,ibox,3,1,hm1)
         call twh_hmatrix(GLB_GET,ibox,3,2,hm2)
         call twh_hmatrix(GLB_GET,ibox,3,3,hm3)
         write(inpnum,'(3f20.12)') hm1,hm2,hm3
c        --- output coordinates
         count_atom = 0
         do imolty = 1, nmolty
            call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
            call twh_nunit(GLB_GET,imolty,nunit)
            do icmt = 1, ncmt
               call twh_chainlist(GLB_GET,icmt,ibox,imolty,ichain)
               do iunit = 1,nunit
                  count_atom = count_atom + 1
                  call twh_ntype(GLB_GET,imolty,iunit,type)
                  call twh_coordinates(GLB_GET,CRD_REAL,ichain
     &                 ,iunit,xtemp,ytemp,ztemp)
c                 --- default format
                  write(inpnum,'(a8,i10)')
     &                 twh_get_nbname(type),count_atom
                  write(inpnum,'(3f20.12)')
     &                 xtemp,ytemp,ztemp
               enddo
            enddo
         enddo
      enddo
c     --- end of CONFIG file

      return
      end
