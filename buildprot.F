#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2011 Marcus G. Martin                       *
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
      subroutine twh_buildprot(lfailure,imolty)
c     ******************************************************************
c     * builds the input file for an amino acid sequence starting from *
c     * the N-terminus                                                 *
c     * imolty: molecule type number                                   *
c     *                                                                *
c     * originally written 02-15-2001 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#include "functions.h"
c     --- debugging flag
#define DEBUG_BUILDPROT 0
c     ---- variables passed to/from the subroutine
      logical lfailure
      integer imolty
c     --- local variables
      logical lfound,lautotorsion

      character*1 stereochem,terminus
      dimension stereochem(NUMAX),terminus(NUMAX)
      character*2 pepname
      dimension pepname(NUMAX)
      character*3 shortname
      character*10 atomname
      dimension atomname(NUMAX)
      character*10 forcefield
      character*15 protgeom,idname
      character*40 aminoname

      integer ipeptide,termflag,iatom,iha,icb,nunit
      integer iprev,total,totatoms,itemp,ivib,iimp,itor
      integer impprev,ndisulfide,myaminonum,myatomnum,ibondend
      dimension myaminonum(NUMAX),myatomnum(NUMAX)
      integer meatom,thematom,me,them,ipair,partner,iother
      integer bondpartner
      dimension bondpartner(NUMAX)

#if DEBUG_BUILDPROT
      write(6,*) 'start BUILDPROT'
#endif
      idname = 'buildprot'

      write(6,*) 'Building the protein input file'
      call twh_read_labeled_string(lfailure,4, 'forcefield', .false.
     &     ,idname,forcefield)
      if ( lfailure ) return

c     --- check for valid forcefield
      if ( forcefield .eq. 'Charmm19' ) then
         write(6,*) '   using the Charmm19 forcefield'
         lautotorsion = .false.
      elseif ( forcefield .eq. 'C19eef1' ) then
         write(6,*) '   using the Charmm19E forcefield'
         lautotorsion = .false.
      elseif ( forcefield .eq. 'Charmm22' ) then
         write(6,*) '   using the Charmm22 forcefield'
         lautotorsion = .true.
      elseif ( forcefield .eq. 'Charmm27' ) then
         write(6,*) '   using the Charmm27 forcefield'
         lautotorsion = .true.
      elseif (forcefield .eq. 'Amber96' ) then
         write(6,*) '   using then Amber-96 forcefield'
         lautotorsion = .true.
      elseif ( forcefield .eq. 'C19sasa' ) then
         write(6,*) '   using the Charmm19S forcefield'
         lautotorsion = .false.
      elseif ( forcefield .eq. 'OPLS-aa' ) then
         write(6,*) '   using the OPLS-aa forcefield'
         lautotorsion = .true.
      else
         write(6,*) 'Invalid forcefield of ',forcefield
         lfailure = .true.
         return
      endif

c     --- read in the protein geometry
      call twh_read_labeled_string(lfailure,4, 'protgeom', .false.
     &     ,idname,protgeom)
      if ( lfailure ) return

c     --- make sure protgeom is valid
      if ( protgeom .eq. 'cyclic' ) then
         write(6,*) '   for a cyclic polypeptide'
      elseif ( protgeom .eq. 'linear' ) then
         write(6,*) '   for a single zwitterionic polypeptide'
         write(6,*) '   or for multiple disulfide bonded zwitterionic'
     &        ,' polypeptides'
      else
         write(6,*) 'protgeom value of ',protgeom,'is not valid'
         lfailure = .true.
         return
      endif

c     --- read in the peptide sequence
      call twh_read_label_formatted(lfailure,4, 'pepname', .false.
     &     , idname, '(a7)')
      if ( lfailure ) return
      call twh_nunit(GLB_GET,imolty,nunit)
      do ipeptide = 1,nunit
         if ( protgeom .eq. 'cyclic' ) then
            read(4,*) pepname(ipeptide),stereochem(ipeptide)
     &           ,bondpartner(ipeptide)
            write(6,15) 'peptide:',ipeptide,'name:',pepname(ipeptide)
     &           ,'stereochem:',stereochem(ipeptide),'bondpartner:'
     &           ,bondpartner(ipeptide)
 15         format(a8,i5,1x,a5,a2,1x,a11,a1,1x,a12,i5)
         elseif ( protgeom .eq. 'linear' ) then
            read(4,*) pepname(ipeptide),stereochem(ipeptide)
     &           ,bondpartner(ipeptide),terminus(ipeptide)
            write(6,16) 'peptide:',ipeptide,'name:',pepname(ipeptide)
     &           ,'stereochem:',stereochem(ipeptide),'bondpartner:'
     &           ,bondpartner(ipeptide),'terminus:',terminus(ipeptide)
 16         format(a8,i5,1x,a5,a2,1x,a11,a1,1x,a12,i5,1x,a9,a1)
         else
            write(6,*) 'BUILDPROT: cannot read sequence for ',protgeom
         endif
      enddo

c     --- initialize all of the bending and torsion arrays
      do iatom = 1,NUMAX
c        --- improper types were already set to 0 in initialize
c        --- intialize disulfide information
         myaminonum(iatom) = 0
         myatomnum(iatom) = 0
      enddo

c     --- initialize iprev and total
      iprev = 1
      total = 0
      totatoms = 0
c     --- initialize ndisulfide
      ndisulfide = 0
c     --- initialize ibondend
      ibondend = 0
c     --- initialize termflag
      termflag = 1

      do ipeptide = 1,nunit
c        --- increment the current atom count by the previous atom count
c        --- plus the number of atoms in the previous residue
         iatom = iprev + total
         if ( protgeom .eq. 'cyclic' ) then
            termflag = 0
         elseif ( protgeom .eq. 'linear' ) then
c           --- determine the terminus from the terminus array
            if ( terminus(ipeptide) .eq. 'N' ) then
c              --- n-term
c              --- means we either need to be at the beginning, or 
c              --- the previous residue was a c-term
               if ( termflag .ne. 1 ) then
                  write(6,*) 'BUILDPROT: problem with terminus'
     &                 ,' for peptide ',ipeptide
                  write(6,*) 'Cannot have an N-terminus unless this'
     &                 ,' is the first peptide'
                  write(6,*) 'or immediately follows a C-term'
                  lfailure = .true.
                  return
               endif
               termflag = -1
            elseif ( terminus(ipeptide) .eq. '-' ) then
c              --- central peptide
c              --- can follow anything except the c-term
               if ( termflag .eq. 1 ) then
                  write(6,*) 'BUILDPROT: problem with terminus'
     &                 ,' for peptide ',ipeptide
                  write(6,*) 'A central peptide cannot immediately'
     &                 ,' follow the c-terminus'
                  write(6,*) 'and it cannot be the first peptide in'
     &                 ,' a chain (use cyclic for that feature)'
                  lfailure = .true.
                  return
               endif
               termflag = 0
            elseif ( terminus(ipeptide) .eq. 'C' ) then
c              --- c-term 
c              --- can follow anything except another c-term
               if ( termflag .eq. 1 ) then
                  write(6,*) 'BUILDPROT: problem with terminus'
     &                 ,' for peptide ',ipeptide
                  write(6,*) 'A c-terminus cannot immediately'
     &                 ,' follow another c-terminus'
                  write(6,*) 'and it cannot be the first peptide in'
     &                 ,'a chain (must list starting from the n-term)'
                  lfailure = .true.
                  return
               endif
               termflag = 1
            else
c              --- unknown termflag
               write(6,*) 'BUILDPROT: unknown termflag'
               write(6,*) terminus(ipeptide)
               lfailure = .true.
               return
            endif
         endif

         if ( forcefield .eq. 'Amber96' ) then
            call twh_polyamber96(lfailure,imolty,iatom,termflag,total
     &           ,ibondend,ipeptide,ndisulfide,myaminonum,myatomnum,iha
     &           ,icb,pepname,atomname,aminoname,shortname)
            if ( lfailure ) return
         elseif ( forcefield .eq. 'Charmm19' ) then
            call twh_polycharmm19(lfailure,imolty,iatom,termflag,total
     &           ,ibondend,ipeptide,iha
     &           ,icb,pepname,atomname,aminoname,shortname)
            if ( lfailure ) return
         elseif ( forcefield .eq. 'C19eef1' ) then
            call twh_polyc19eef1(lfailure,imolty,iatom,termflag,total
     &           ,ibondend,ipeptide,iha
     &           ,icb,pepname,atomname,aminoname,shortname)
            if ( lfailure ) return
         elseif ( forcefield .eq. 'Charmm22' ) then
            call twh_polycharmm22(lfailure,imolty,iatom,termflag,total
     &           ,ibondend,ipeptide,ndisulfide,myaminonum,myatomnum,iha
     &           ,icb,pepname,atomname,aminoname,shortname)
            if ( lfailure ) return
         elseif ( forcefield .eq. 'Charmm27' ) then
            call twh_polycharmm27(lfailure,imolty,iatom,termflag,total
     &           ,ibondend,ipeptide,ndisulfide,myaminonum,myatomnum,iha
     &           ,icb,pepname,atomname,aminoname,shortname)
            if ( lfailure ) return
         elseif (forcefield .eq. 'C19sasa') then
            call twh_polyc19eef1(lfailure,imolty,iatom,termflag,total
     &           ,ibondend,ipeptide,iha
     &           ,icb,pepname,atomname,aminoname,shortname)
            if ( lfailure ) return
         elseif (forcefield .eq. 'OPLS-aa') then
            call twh_polyoplsaa(lfailure,imolty,iatom,termflag,total
     &           ,ibondend,ipeptide,ndisulfide,myaminonum,myatomnum,iha
     &           ,icb,pepname,atomname,aminoname,shortname)
            if ( lfailure ) return
         else
            write(6,*) 'BUILDPROT: unknown forcefield'
            write(6,*) forcefield
            lfailure = .true.
            return
         endif

c        --- set up the impropers for the stereochemistry
c        --- check stereochem as only racemic is valid for glycine
         if ( shortname .eq. 'GLY' .and.
     &        .not.(stereochem(ipeptide) .eq. 'r'
     &        .or. stereochem(ipeptide) .eq. 'R') ) then
            write(6,*) 'stereochem of glycine must be (r) racemic'
            write(6,*) 'please change invalid stereochem of '
     &           ,stereochem(ipeptide)
            lfailure = .true.
            return
         elseif ( shortname .eq. 'FOR' .and.
     &        .not.(stereochem(ipeptide) .eq. 'r'
     &        .or. stereochem(ipeptide) .eq. 'R') ) then
            write(6,*) 'stereochem of formaldehyde cap must be (r)'
     &           ,' racemic'
            write(6,*) 'please change invalid stereochem of '
     &           ,stereochem(ipeptide)
            lfailure = .true.
            return
         elseif ( shortname .eq. 'ACE' .and.
     &        .not.(stereochem(ipeptide) .eq. 'r'
     &        .or. stereochem(ipeptide) .eq. 'R') ) then
            write(6,*) 'stereochem of ACE cap must be (r)'
     &           ,' racemic'
            write(6,*) 'please change invalid stereochem of '
     &           ,stereochem(ipeptide)
            lfailure = .true.
            return
         elseif ( shortname .eq. 'NHE' .and.
     &        .not.(stereochem(ipeptide) .eq. 'r'
     &        .or. stereochem(ipeptide) .eq. 'R') ) then
            write(6,*) 'stereochem of NHE cap must be (r)'
     &           ,' racemic'
            write(6,*) 'please change invalid stereochem of '
     &           ,stereochem(ipeptide)
            lfailure = .true.
            return
         endif

         if ( forcefield .eq. 'Charmm22'
     &        .or. forcefield .eq. 'Amber96'
     &        .or. forcefield .eq. 'Charmm27') then
            if ( stereochem(ipeptide) .eq. 'l'
     &           .or. stereochem(ipeptide) .eq. 'L' ) then
               write(6,*) 'L-',aminoname
               impprev = twh_get_inimprop(imolty,iatom+1)
               call twh_inimprop(GLB_SET,imolty,iatom+1,3+impprev)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,impprev+1
     &              ,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,impprev+1,iha)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,impprev+1,icb)
               call twh_itimprop(GLB_SET,imolty,iatom+1,impprev+1,-1)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,impprev+2
     &              ,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,impprev+2
     &              ,iatom+2)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,impprev+2,iha)
               call twh_itimprop(GLB_SET,imolty,iatom+1,impprev+2,-1)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,impprev+3
     &              ,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,impprev+3,icb)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,impprev+3
     &              ,iatom+2)
               call twh_itimprop(GLB_SET,imolty,iatom+1,impprev+3,-1)
            elseif ( stereochem(ipeptide) .eq. 'd'
     &              .or. stereochem(ipeptide) .eq. 'D' ) then
               write(6,*) 'D-',aminoname
               impprev = twh_get_inimprop(imolty,iatom+1)
               call twh_inimprop(GLB_SET,imolty,iatom+1,3+impprev)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,impprev+1
     &              ,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,impprev+1,icb)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,impprev+1,iha)
               call twh_itimprop(GLB_SET,imolty,iatom+1,impprev+1,-1)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,impprev+2
     &              ,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,impprev+2,iha)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,impprev+2
     &              ,iatom+2)
               call twh_itimprop(GLB_SET,imolty,iatom+1,impprev+2,-1)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,impprev+3
     &              ,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,impprev+3
     &              ,iatom+2)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,impprev+3,icb)
               call twh_itimprop(GLB_SET,imolty,iatom+1,impprev+3,-1)
            elseif ( stereochem(ipeptide) .eq. 'r'
     &              .or. stereochem(ipeptide) .eq. 'R' ) then
               write(6,*) 'racemic ',aminoname
            else
               write(6,*) 'stereochem of ',stereochem(ipeptide)
     &              ,' is not valid'
               lfailure = .true.
               return
            endif
         else
            write(6,*) aminoname
         endif
         totatoms = totatoms + total

         if ( totatoms .gt. NUMAX ) then
            write(6,*) 'total, totatoms',total,totatoms
            write(6,*) 'totatoms > NUMAX'
            lfailure = .true.
            return
         endif

c        --- assign the amino acid shortcut names to the atoms
         do itemp = iatom,iatom+total-1
            call twh_aminoshort(GLB_SET,imolty,itemp,shortname)
            call twh_aminonum(GLB_SET,imolty,itemp,ipeptide)
         enddo

c        --- update iprev
         iprev = iatom
         if ( iprev +35 .gt. NUMAX ) then
            write(6,*) 'BUILDPROT: atom count approaching NUMAX'
            write(6,*) 'current atom count',totatoms
            write(6,*) 'peptide count:',ipeptide
         endif
      enddo

c     --- fix up problems if this is a cyclic polypeptide
      if ( protgeom .eq. 'cyclic' ) then
c        --- search through the vibration and improper list
c        --- reset any bonds with atom 0 to iatom+2
c        --- reset any bonds with iatom+total to 1
         do itemp = 1,totatoms
            do ivib = 1,twh_get_invib(imolty,itemp)
               if ( twh_get_ijvib(imolty,itemp,ivib) .eq. iatom+total 
     &              ) then
                  call twh_ijvib(GLB_SET,imolty,itemp,ivib,1)
               elseif ( twh_get_ijvib(imolty,itemp,ivib) .eq. 0 ) then
                  call twh_ijvib(GLB_SET,imolty,itemp,ivib,iatom+2)
               endif
            enddo
            if (forcefield .eq. 'Charmm19' .or.
     &          forcefield .eq. 'C19eef1' .or. 
     &          forcefield .eq. 'C19sasa') then
               do itor = 1,twh_get_intor(imolty,itemp)
                  if ( twh_get_ijtor2(imolty,itemp,itor) 
     &                 .eq. iatom+total)
     &                 call twh_ijtor2(GLB_SET,imolty,itemp,itor,1)
                  if ( twh_get_ijtor3(imolty,itemp,itor) 
     &                 .eq. iatom+total)
     &                 call twh_ijtor3(GLB_SET,imolty,itemp,itor,1)
                  if ( twh_get_ijtor4(imolty,itemp,itor) 
     &                 .eq. iatom+total)
     &                 call twh_ijtor4(GLB_SET,imolty,itemp,itor,1)
                  if ( twh_get_ijtor2(imolty,itemp,itor) .eq. 0)
     &                 call twh_ijtor2(GLB_SET,imolty,itemp,itor
     &                 ,iatom+2)
                  if ( twh_get_ijtor3(imolty,itemp,itor) .eq. 0)
     &                 call twh_ijtor3(GLB_SET,imolty,itemp,itor
     &                 ,iatom+2)
                  if ( twh_get_ijtor4(imolty,itemp,itor) .eq. 0)
     &                 call twh_ijtor4(GLB_SET,imolty,itemp,itor
     &                 ,iatom+2)
                  if ( twh_get_ijtor2(imolty,itemp,itor) 
     &                 .eq. iatom+total+1)
     &                 call twh_ijtor2(GLB_SET,imolty,itemp,itor,2)
                  if ( twh_get_ijtor3(imolty,itemp,itor) 
     &                 .eq. iatom+total+1)
     &                 call twh_ijtor3(GLB_SET,imolty,itemp,itor,2)
                  if ( twh_get_ijtor4(imolty,itemp,itor) 
     &                 .eq. iatom+total+1)
     &                 call twh_ijtor4(GLB_SET,imolty,itemp,itor,2)
                  if ( twh_get_ijtor2(imolty,itemp,itor) .eq. -1)
     &                 call twh_ijtor2(GLB_SET,imolty,itemp,itor
     &                 ,iatom+1)
                  if ( twh_get_ijtor3(imolty,itemp,itor) .eq. -1)
     &                 call twh_ijtor3(GLB_SET,imolty,itemp,itor
     &                 ,iatom+1)
                  if ( twh_get_ijtor4(imolty,itemp,itor) .eq. -1)
     &                 call twh_ijtor4(GLB_SET,imolty,itemp,itor
     &                 ,iatom+1)
                  if ( twh_get_ijtor2(imolty,itemp,itor) 
     &                 .eq. iatom+total+2)
     &                 call twh_ijtor2(GLB_SET,imolty,itemp,itor,3)
                  if ( twh_get_ijtor3(imolty,itemp,itor) 
     &                 .eq. iatom+total+2)
     &                 call twh_ijtor3(GLB_SET,imolty,itemp,itor,3)
                  if ( twh_get_ijtor4(imolty,itemp,itor) 
     &                 .eq. iatom+total+2)
     &                 call twh_ijtor4(GLB_SET,imolty,itemp,itor,3)
                  if ( twh_get_ijtor2(imolty,itemp,itor) .eq. -2)
     &                 call twh_ijtor2(GLB_SET,imolty,itemp,itor,iatom)
                  if ( twh_get_ijtor3(imolty,itemp,itor) .eq. -2)
     &                 call twh_ijtor3(GLB_SET,imolty,itemp,itor,iatom)
                  if ( twh_get_ijtor4(imolty,itemp,itor) .eq. -2)
     &                 call twh_ijtor4(GLB_SET,imolty,itemp,itor,iatom)
               enddo
            endif
            do iimp = 1,twh_get_inimprop(imolty,itemp)
               if ( twh_get_ijimprop2(imolty,itemp,iimp) 
     &              .eq. iatom+total)
     &              call twh_ijimprop2(GLB_SET,imolty,itemp,iimp,1)
               if ( twh_get_ijimprop3(imolty,itemp,iimp) 
     &              .eq. iatom+total)
     &              call twh_ijimprop3(GLB_SET,imolty,itemp,iimp,1)
               if ( twh_get_ijimprop4(imolty,itemp,iimp) 
     &              .eq. iatom+total)
     &              call twh_ijimprop4(GLB_SET,imolty,itemp,iimp,1)
               if ( twh_get_ijimprop2(imolty,itemp,iimp) .eq. 0)
     &              call twh_ijimprop2(GLB_SET,imolty,itemp,iimp
     &              ,iatom+2)
               if ( twh_get_ijimprop3(imolty,itemp,iimp) .eq. 0)
     &              call twh_ijimprop3(GLB_SET,imolty,itemp,iimp
     &              ,iatom+2)
               if ( twh_get_ijimprop4(imolty,itemp,iimp) .eq. 0)
     &              call twh_ijimprop4(GLB_SET,imolty,itemp,iimp
     &              ,iatom+2)
            enddo
         enddo
      endif

c     --- attach all of the disulfide bonds to each other
      ipair = 0
      do while ( ndisulfide .ne. 0 )
         me = myaminonum(1)
         partner = bondpartner(me)
         lfound = .false.
         do iother = 2,ndisulfide
            if ( myaminonum(iother) .eq. partner ) then
               lfound = .true.
               ipair = iother
               them = myaminonum(iother)
c              --- this is my partner, make sure reverse is true
               if ( bondpartner(them) .ne. me ) then
                  write(6,*) '1 peptide ',me,' bondpartner ',partner
                  write(6,*) '2 peptide ',them
     &                 ,' bondpartner ',bondpartner(them)
                  write(6,*) 'should have matched'
                  lfailure = .true.
                  return
               endif
            endif
         enddo
         if ( .not. lfound ) then
            write(6,*) 'peptide ',me,' bondpartner ',partner
     &           ,' was not found '
            lfailure = .true.
            return
         endif
c        --- hook them together
         meatom = myatomnum(1)
         thematom = myatomnum(ipair)
         call twh_ijvib(GLB_SET,imolty,meatom,2,thematom)
         call twh_ijvib(GLB_SET,imolty,thematom,2,meatom)
c        --- remove the atoms from the nsulfide data structure
         do iother = 1,ndisulfide
            if ( iother+1 .ge. ipair ) then
               myaminonum(iother) = myaminonum(iother+2)
               myatomnum(iother) = myatomnum(iother+2)
            else
               myaminonum(iother) = myaminonum(iother+1)
               myatomnum(iother) = myatomnum(iother+1)
            endif
         enddo
c        --- decrement ndisulfide
         ndisulfide = ndisulfide - 2
      enddo
      call twh_assemble(lfailure,imolty,totatoms,0,forcefield,atomname
     &     ,'standard','manual',lautotorsion)
      if ( lfailure ) return
#if DEBUG_BUILDPROT
      write(6,*) 'finished BUILDPROT'
#endif
      return
      end
