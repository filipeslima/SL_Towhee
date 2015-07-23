#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000-2014 Marcus G. Martin                       *
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
      subroutine twh_checkstruc(lfailure)
c     ******************************************************************
c     * this subroutine runs through the input file structure and makes*
c     * sure that all of the bonds, bends, and torsions are set up in a*
c     * consistent way.                                                *
c     * it then sets up some arrays of information that is useful when *
c     * computing bending angles, and torsions                         *
c     *                                                                *
c     * Originally written 12-15-2000 by M.G. Martin                   *
c     * Last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_ANGLESTYLE
#define FUNCTION_GET_BENCOEFF
#define FUNCTION_GET_BONDPATT
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#define FUNCTION_GET_IBTBEN1
#define FUNCTION_GET_IBTBEN2
#define FUNCTION_GET_IJAA0
#define FUNCTION_GET_IJAA1
#define FUNCTION_GET_IJAA2
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#define FUNCTION_GET_IJOF5
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INBNBN
#define FUNCTION_GET_INOF
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_ITAA
#define FUNCTION_GET_ITBEN
#define FUNCTION_GET_ITVIB
#define FUNCTION_GET_INAA
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_ITIMPROP
#define FUNCTION_GET_ITTOR
#define FUNCTION_GET_LBONBON
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_GET_VIBCOEFF
#include "functions.h"
#define DEBUG_CHECKSTRUC 0
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- local logicals
      logical lfound,lbnbn,lintraterm,loutdft,lonefourtor
c     --- integer scalars
      integer imolty,ivib,jvib,iben,jben,itor,jtor,iunit,iaa,jaa,itype
      integer iimp,jimp,itest,itestb,itestc,junit,count,iof,jof,itemp
      integer aunit,bunit,cunit,kunit,atype,ii,kk,ktype,matom,type
      integer iorder,nmolty,nunit,nstep,numtor
c     --- local double precision arrays
      integer torcross
      dimension torcross(NUMAX)
c     --- local dp scalars
      double precision rij,rjk,rik,zi,zk,theta,dvalue

      write(6,*) 'Verifying input structures are consistent'

c     --- initialize flag to tell if we have any bond-bond interactions
      lbnbn = .false.
c     --- determine whether we need to stop if we observe a non-valid
c     --- intrmolecular setup.  we might not want to stop if all we are
c     --- doing is converting from lammps into tramonto and not doing
c     --- any towhee simulation steps
      call twh_nstep(GLB_GET,nstep)
      call twh_loutdft(GLB_GET,loutdft)
      if ( loutdft .and. nstep .eq. 0 ) then
         lintraterm = .false.
      else
         lintraterm = .true.
      endif

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1, nmolty
c        --- check that each atom has self-consistent and valid values
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1, nunit
c           --- check vibrations
            do ivib = 1, twh_get_invib(imolty,iunit)
c              --- check that vibration atoms are between 1 and nunit
               if ( twh_get_ijvib(imolty,iunit,ivib) .lt. 1 .or.
     &              twh_get_ijvib(imolty,iunit,ivib) .gt. nunit
     &              )then
                  write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                 ,iunit,ivib
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,twh_get_ijvib(imolty,iunit,ivib)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that vibration type is between 1 and tvibmax
               if ( twh_get_itvib(imolty,iunit,ivib) .lt. 1 .or.
     &              twh_get_itvib(imolty,iunit,ivib) .gt. TVIBMAX )then
                  write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                 ,iunit,ivib
                  write(6,*) 'Vibration type is nonvalid value '
     &                 ,twh_get_itvib(imolty,iunit,ivib)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
c           --- check for multiple vibrations across the same atoms
            do ivib = 1, twh_get_invib(imolty,iunit)-1
               do jvib = ivib+1, twh_get_invib(imolty,iunit)
                  if ( twh_get_ijvib(imolty,iunit,ivib) .eq. 
     &                 twh_get_ijvib(imolty,iunit,jvib) ) then
                     write(6,*) 'Problem: imolty,iunit,ivib,jvib',
     &                    imolty,iunit,ivib,jvib
                     write(6,*) 'Are both connected to atom '
     &                    ,twh_get_ijvib(imolty,iunit,ivib)
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo
c           --- check that each pair of atoms that is bonded to this
c           --- atom also have a bending interaction
            do ivib = 1, twh_get_invib(imolty,iunit)-1
               aunit = twh_get_ijvib(imolty,iunit,ivib)
               do jvib = ivib+1,twh_get_invib(imolty,iunit)
                  bunit = twh_get_ijvib(imolty,iunit,jvib)
                  lfound = .false.
                  do iben = 1,twh_get_inben(imolty,aunit)
                     if ( twh_get_ijben2(imolty,aunit,iben) .eq. iunit 
     &                    .and.
     &                    twh_get_ijben3(imolty,aunit,iben) .eq. bunit 
     &                    ) then
                        lfound = .true.
                     endif
                  enddo
                  if ( .not. lfound ) then
                     write(6,*) 'imolty,iunit,ivib,jvib'
     &                    ,imolty,iunit,ivib,jvib
                     write(6,*) 'atoms ',aunit,bunit,' are connected to'
     &                    ,' iunit, but have no bend'
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo

c           --- check bending
            do iben = 1, twh_get_inben(imolty,iunit)
c              --- check that bending atoms are between 1 and nunit
               if ( twh_get_ijben2(imolty,iunit,iben) .lt. 1 .or.
     &              twh_get_ijben2(imolty,iunit,iben) .gt. nunit
     &              .or.
     &              twh_get_ijben3(imolty,iunit,iben) .lt. 1 .or.
     &              twh_get_ijben3(imolty,iunit,iben) .gt. nunit
     &              )then
                  write(6,*) 'Problem: imolty,iunit,iben',imolty
     &                 ,iunit,iben
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,twh_get_ijben2(imolty,iunit,iben)
     &                 ,twh_get_ijben3(imolty,iunit,iben)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that bending type is between 1 and tbenmax
               if ( twh_get_itben(imolty,iunit,iben) .lt. 1 .or.
     &              twh_get_itben(imolty,iunit,iben) .gt. TBENMAX 
     &              )then
                  write(6,*) 'Problem: imolty,iunit,iben',imolty
     &                 ,iunit,iben
                  write(6,*) 'Bending type is nonvalid value '
     &                 ,twh_get_itben(imolty,iunit,iben)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
c           --- check for multiple bends across the same atoms
            do iben = 1, twh_get_inben(imolty,iunit)-1
               do jben = iben+1, twh_get_inben(imolty,iunit)
                  if ( twh_get_ijben2(imolty,iunit,iben) .eq. 
     &                 twh_get_ijben2(imolty,iunit,jben) .and. 
     &                 twh_get_ijben3(imolty,iunit,iben) .eq. 
     &                 twh_get_ijben3(imolty,iunit,jben) ) then
                     write(6,*) 'Problem: imolty,iunit,iben,jben',
     &                    imolty,iunit,iben,jben
                     write(6,*) 'Are both connected to atoms '
     &                    ,twh_get_ijben2(imolty,iunit,iben)
     &                    ,twh_get_ijben3(imolty,iunit,iben)
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo
c           --- check that all bends which are also bonded to another
c           --- atom have a torsion interaction.
            do iben = 1, twh_get_inben(imolty,iunit)
               aunit = twh_get_ijben3(imolty,iunit,iben)
               bunit = twh_get_ijben2(imolty,iunit,iben)
               do ivib = 1, twh_get_invib(imolty,iunit)
                  cunit = twh_get_ijvib(imolty,iunit,ivib)
                  if ( cunit .ne. bunit .and. cunit .ne. aunit ) then
                     lfound = .false.
                     call twh_intor(GLB_GET,imolty,cunit,numtor)
                     do itor = 1,numtor
                        if ( twh_get_ijtor2(imolty,cunit,itor).eq.iunit
     &                       .and.
     &                       twh_get_ijtor3(imolty,cunit,itor).eq.bunit
     &                       .and.
     &                       twh_get_ijtor4(imolty,cunit,itor).eq.aunit
     &                       ) then
                           lfound = .true.
                        endif
                     enddo
                     if ( .not. lfound ) then
                        write(6,*) 'imolty,iunit,iben,ivib',imolty,iunit
     &                       ,iben,ivib
                        write(6,*) 'bonded to an atom and an angle that'
     &                       ,' do not form a torsion'
                        if ( lintraterm ) then
                           lfailure = .true.
                           return
                        endif
                     endif
                  endif
               enddo
            enddo

c           --- check torsions
            call twh_intor(GLB_GET,imolty,iunit,numtor)
            do itor = 1, numtor
c              --- check that torsion atoms are between 1 and nunit
               if ( twh_get_ijtor2(imolty,iunit,itor) .lt. 1 .or.
     &              twh_get_ijtor2(imolty,iunit,itor) .gt. nunit
     &              .or.
     &              twh_get_ijtor3(imolty,iunit,itor) .lt. 1 .or.
     &              twh_get_ijtor3(imolty,iunit,itor) .gt. nunit
     &              .or.
     &              twh_get_ijtor4(imolty,iunit,itor) .lt. 1 .or.
     &              twh_get_ijtor4(imolty,iunit,itor) .gt. nunit
     &              )then
                  write(6,*) 'Problem: imolty,iunit,itor',imolty
     &                 ,iunit,itor
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,twh_get_ijtor2(imolty,iunit,itor)
     &                 ,twh_get_ijtor3(imolty,iunit,itor)
     &                 ,twh_get_ijtor4(imolty,iunit,itor)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif

c              --- check that torsion type is between 1 and ttormax
               if ( abs(twh_get_ittor(imolty,iunit,itor)) .eq. 0 .or.
     &              abs(twh_get_ittor(imolty,iunit,itor)) .gt. TTORMAX 
     &              )then
                  write(6,*) 'Problem: imolty,iunit,itor',imolty
     &                 ,iunit,itor
                  write(6,*) 'Torsion type is nonvalid value '
     &                 ,twh_get_ittor(imolty,iunit,itor)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo

c           --- check for multiple torsions across the same atoms
            call twh_intor(GLB_GET,imolty,iunit,numtor)
            do itor = 1, numtor-1
               do jtor = itor+1, twh_get_intor(imolty,iunit)
                  if ( twh_get_ijtor2(imolty,iunit,itor) .eq. 
     &                 twh_get_ijtor2(imolty,iunit,jtor) .and. 
     &                 twh_get_ijtor3(imolty,iunit,itor) .eq. 
     &                 twh_get_ijtor3(imolty,iunit,jtor) .and.
     &                 twh_get_ijtor4(imolty,iunit,itor) .eq. 
     &                 twh_get_ijtor4(imolty,iunit,jtor) ) then
                     write(6,*) 'Problem: imolty,iunit,itor,jtor',
     &                    imolty,iunit,itor,jtor
                     write(6,*) 'Are both connected to atoms '
     &                    ,twh_get_ijtor2(imolty,iunit,itor)
     &                    ,twh_get_ijtor3(imolty,iunit,itor)
     &                    ,twh_get_ijtor4(imolty,iunit,itor)
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo

c           --- check angle-angle terms
            do iaa = 1, twh_get_inaa(imolty,iunit)
c              --- check that ang-ang atoms are between 1 and nunit
               if ( twh_get_ijaa0(imolty,iunit,iaa) .lt. 1 .or.
     &              twh_get_ijaa0(imolty,iunit,iaa) .gt. nunit
     &              .or.
     &              twh_get_ijaa1(imolty,iunit,iaa) .lt. 1 .or.
     &              twh_get_ijaa1(imolty,iunit,iaa) .gt. nunit
     &              .or.
     &              twh_get_ijaa2(imolty,iunit,iaa) .lt. 1 .or.
     &              twh_get_ijaa2(imolty,iunit,iaa) .gt. nunit
     &              )then
                  write(6,*) 'Problem: imolty,iunit,iaa',imolty
     &                 ,iunit,iaa
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,twh_get_ijaa0(imolty,iunit,iaa)
     &                 ,twh_get_ijaa1(imolty,iunit,iaa)
     &                 ,twh_get_ijaa2(imolty,iunit,iaa)
                  if ( lintraterm )  then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that ang-ang type is between 1 and taamax
               if ( twh_get_itaa(imolty,iunit,iaa) .lt. 1 .or.
     &              twh_get_itaa(imolty,iunit,iaa) .gt. TAAMAX )then
                  write(6,*) 'Problem: imolty,iunit,iaa',imolty
     &                 ,iunit,iaa
                  write(6,*) 'Ang-ang type is nonvalid value '
     &                 ,twh_get_itaa(imolty,iunit,iaa)
                  if ( lintraterm )  then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
c           --- check for multiple ang-ang across the same atoms
            do iaa = 1, twh_get_inaa(imolty,iunit)-1
               do jaa = iaa+1, twh_get_inaa(imolty,iunit)
                  if ( twh_get_ijaa0(imolty,iunit,iaa) .eq. 
     &                 twh_get_ijaa0(imolty,iunit,jaa) .and. 
     &                 twh_get_ijaa1(imolty,iunit,iaa) .eq. 
     &                 twh_get_ijaa1(imolty,iunit,jaa) .and.
     &                 twh_get_ijaa2(imolty,iunit,iaa) .eq. 
     &                 twh_get_ijaa2(imolty,iunit,jaa) ) then
                     write(6,*) 'Problem: imolty,iunit,iaa,jaa',
     &                    imolty,iunit,iaa,jaa
                     write(6,*) 'Are both connected to atoms '
     &                    ,twh_get_ijaa0(imolty,iunit,iaa)
     &                    ,twh_get_ijaa1(imolty,iunit,iaa)
     &                    ,twh_get_ijaa2(imolty,iunit,iaa)
                     if ( lintraterm )  then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo

c           --- check improper torsions
            do iimp = 1, twh_get_inimprop(imolty,iunit)
c              --- check that atoms are between 1 and nunit
               if ( twh_get_ijimprop2(imolty,iunit,iimp) .lt. 1 .or.
     &              twh_get_ijimprop2(imolty,iunit,iimp) 
     &              .gt. nunit .or.
     &              twh_get_ijimprop3(imolty,iunit,iimp) .lt. 1 .or.
     &              twh_get_ijimprop3(imolty,iunit,iimp) 
     &              .gt. nunit .or.
     &              twh_get_ijimprop4(imolty,iunit,iimp) .lt. 1 .or.
     &              twh_get_ijimprop4(imolty,iunit,iimp) 
     &              .gt. nunit)then
                  write(6,*) 'Problem: imolty,iunit,iimp',imolty
     &                 ,iunit,iimp
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,twh_get_ijimprop2(imolty,iunit,iimp)
     &                 ,twh_get_ijimprop3(imolty,iunit,iimp)
     &                 ,twh_get_ijimprop4(imolty,iunit,iimp)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that imp torsion type is between 1 and timpmax
               if ( twh_get_itimprop(imolty,iunit,iimp) .lt. 1 .or.
     &              twh_get_itimprop(imolty,iunit,iimp) .gt. TIMPMAX 
     &              )then
                  write(6,*) 'Problem: imolty,iunit,iimp',imolty
     &                 ,iunit,iimp
                  write(6,*) 'ijimprop2,ijimprop3,ijimprop4'
     &                 ,twh_get_ijimprop2(imolty,iunit,iimp)
     &                 ,twh_get_ijimprop3(imolty,iunit,iimp)
     &                 ,twh_get_ijimprop4(imolty,iunit,iimp)
                  write(6,*) 'Improper torsion type is nonvalid value '
     &                 ,twh_get_itimprop(imolty,iunit,iimp)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
c           --- check for multiple imp torsions across the same atoms
            do iimp = 1, twh_get_inimprop(imolty,iunit)-1
               do jimp = iimp+1, twh_get_inimprop(imolty,iunit)
                  if ( twh_get_ijimprop2(imolty,iunit,iimp) .eq. 
     &                 twh_get_ijimprop2(imolty,iunit,jimp) .and. 
     &                 twh_get_ijimprop3(imolty,iunit,iimp) .eq. 
     &                 twh_get_ijimprop3(imolty,iunit,jimp) .and.
     &                 twh_get_ijimprop4(imolty,iunit,iimp) .eq. 
     &                 twh_get_ijimprop4(imolty,iunit,jimp) ) then
                     write(6,*) 'Problem: imolty,iunit,iimp,jimp',
     &                    imolty,iunit,iimp,jimp
                     write(6,*) 'Are both connected to atoms '
     &                    ,twh_get_ijimprop2(imolty,iunit,iimp)
     &                    ,twh_get_ijimprop3(imolty,iunit,iimp)
     &                    ,twh_get_ijimprop4(imolty,iunit,iimp)
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo
c           --- check one-five interactions
            do iof = 1,twh_get_inof(imolty,iunit)
               matom = twh_get_ijof5(imolty,iunit,iof)
               if ( matom .lt. 1 .or. matom .gt. nunit ) then
c                 --- ijof5 out of bounds
                  write(6,*) 'CHECKSTRUC: atom out of allowed range'
                  write(6,*) 'imolty,iunit,ijof5',imolty,iunit,matom
                  lfailure = .true.
                  return
               endif
               call twh_itof(GLB_GET,imolty,iunit,iof,type)
               if ( type .lt. 1 .or. type .gt. TOFMAX ) then
c                 --- itof out of bound
                  write(6,*) 'CHECKSTRUC: type out of allowed range'
                  write(6,*) 'imolty,iunit,itof',imolty,iunit,type
                  lfailure = .true.
                  return
               endif
c              --- check multiples
               do jof = 1,twh_get_inof(imolty,iunit)
                  if ( iof .ne. jof ) then
                     if ( twh_get_ijof5(imolty,iunit,jof) 
     &                    .eq. matom ) then
                        write(6,*) 'CHECKSTRUC: multiple identical 1-5'
                        write(6,*) 'imolty,iunit,iof,ijof5',imolty,iunit
     &                       ,iof,twh_get_ijof5(imolty,iunit,iof)
                        write(6,*) 'imolty,iunit,jof,ijof5',imolty,iunit
     &                       ,jof,twh_get_ijof5(imolty,iunit,jof)
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo
         enddo

c        --- check the terms are consistent from one atom to the next
c        --- for the terms that should be symmetric
         do iunit = 1,nunit
c           --- check vibrations
            do ivib = 1, twh_get_invib(imolty,iunit)
c              --- check that atom is not connected to itself
               if ( twh_get_ijvib(imolty,iunit,ivib) .eq. iunit )then
                  write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                 ,iunit,ivib
                  write(6,*) 'Atom is connected to itself '
     &                 ,twh_get_ijvib(imolty,iunit,ivib)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that each other atom also contains this 
c              --- interaction
               itest = twh_get_ijvib(imolty,iunit,ivib)
               lfound = .false.
               do jvib = 1, twh_get_invib(imolty,itest)
                  if ( twh_get_ijvib(imolty,itest,jvib) .eq. iunit 
     &                 ) then
                     if ( twh_get_itvib(imolty,itest,jvib) .ne. 
     &                    twh_get_itvib(imolty,iunit,ivib) ) then
                        write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                       ,iunit,ivib
                        write(6,*) 'Has nonsymmetric type: itest,jvib'
     &                       ,itest,jvib
                        write(6,*) 'iunit ivib type '
     &                       ,twh_get_itvib(imolty,iunit,ivib)
                        write(6,*) 'itest jvib type '
     &                       ,twh_get_itvib(imolty,itest,jvib)
                        if ( lintraterm ) then
                           lfailure = .true.
                           return
                        endif
                     endif
                     lfound = .true.
                  endif
               enddo
               if ( .not. lfound ) then 
                  write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                 ,iunit,ivib
                  write(6,*) 'There is no symmetric vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
         enddo

         do iunit = 1,nunit
c           --- check bends
            do iben = 1, twh_get_inben(imolty,iunit)
               itest = twh_get_ijben2(imolty,iunit,iben)
               itestb = twh_get_ijben3(imolty,iunit,iben)
c              --- check that atom #2 is bonded to this atom
               lfound = .false.
               do ivib = 1,twh_get_invib(imolty,iunit)
                  if ( twh_get_ijvib(imolty,iunit,ivib) .eq. itest 
     &                 ) then
                     lfound = .true.
c                    --- store bond type in this angle
                     itemp = twh_get_itvib(imolty,iunit,ivib)
                     call twh_ibtben1(GLB_SET,imolty,iunit,iben,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijben2,ijben3,iben'
     &                 ,imolty,iunit,itest,itestb,iben
                  write(6,*) 'Atoms 1-2 connected in a bend are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #3 is bonded to atom #2
               lfound = .false.
               do ivib = 1,twh_get_invib(imolty,itest)
                  if ( twh_get_ijvib(imolty,itest,ivib) .eq. itestb 
     &                 ) then
                     lfound = .true.
c                    --- store bond type in this angle
                     itemp = twh_get_itvib(imolty,itest,ivib)
                     call twh_ibtben2(GLB_SET,imolty,iunit,iben,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijben2,ijben3,iben'
     &                 ,imolty,iunit,itest,itestb,iben
                  write(6,*) 'Atoms 2-3 connected in a bend are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif

c              --- check that each other atom also contains this 
c              --- interaction
               lfound = .false.
               do jben = 1, twh_get_inben(imolty,itestb)
                  if ( twh_get_ijben3(imolty,itestb,jben) .eq. iunit 
     &                 .and. 
     &                 twh_get_ijben2(imolty,itestb,jben) .eq. itest 
     &                 ) then
                     if ( twh_get_itben(imolty,itestb,jben) .ne. 
     &                    twh_get_itben(imolty,iunit,iben) ) then
                        write(6,*) 'Problem: imolty,iunit,iben',imolty
     &                       ,iunit,iben
                        write(6,*) 'Has nonsymmetric type: itestb,itest'
     &                       ,',jben',itestb,itest,jben
                        write(6,*) 'iunit iben type '
     &                       ,twh_get_itben(imolty,iunit,iben)
                        write(6,*) 'itestb jben type '
     &                       ,twh_get_itben(imolty,itestb,jben)
                        if ( lintraterm ) then
                           lfailure = .true.
                           return
                        endif
                     endif
                     lfound = .true.
                  endif
               enddo
               if ( .not. lfound ) then 
                  write(6,*) 'Problem: imolty,iunit,iben',imolty
     &                 ,iunit,iben
                  write(6,*) 'There is no symmetric bend'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif

               if ( twh_get_lbonbon(
     &              twh_get_itben(imolty,iunit,iben)) ) then
c                 --- we have at least on bond-bond interaction
                  lbnbn = .true.
c                 --- do not double count as we are on the center atom
                  if ( iunit .lt. itestb ) then
                     call twh_inbnbn(GLB_INCR,imolty,itest,1)
                     count = twh_get_inbnbn(imolty,itest)
                     itemp = twh_get_itben(imolty,iunit,iben)
                     call twh_orderben(GLB_GET,imolty,iunit,iben,iorder)
                     if ( iorder .eq. ORD_REVERSE ) then
                        itemp = -itemp
                     endif
                     call twh_itbnbn(GLB_SET,imolty,itest,count,itemp)
                     call twh_ijbnbn1(GLB_SET,imolty,itest,count,iunit)
                     call twh_ijbnbn2(GLB_SET,imolty,itest,count,itestb)
                     itemp = twh_get_ibtben1(imolty,iunit,iben)
                     call twh_ibtbnbn1(GLB_SET,imolty,itest,count,itemp)
                     itemp = twh_get_ibtben2(imolty,iunit,iben)
                     call twh_ibtbnbn2(GLB_SET,imolty,itest,count,itemp)
                  endif
               endif
               atype = twh_get_itben(imolty,iunit,iben)
               if ( twh_get_anglestyle(atype) .eq. ANG_STYLE_FOURIER
     &              .or. twh_get_anglestyle(atype)
     &              .eq. ANG_STYLE_TTFOURIER ) then
c                 --- this only works with the UFF 12-6 classical
c                 --- potential
                  if ( twh_get_classical_potential() 
     &                 .ne. 'UFF 12-6' ) then
                     write(6,*) 'CHECKSTRUC: cannot combine anglestyle'
     &                    ,twh_get_anglestyle(atype)
                     write(6,*) 'with anything other than a'
     &                    ,' classical_potential of UFF 12-6'
                     lfailure = .true.
                     return
                  endif
c                 --- using the UFF angle terms means we need to 
c                 --- setup the prefactor here to save time later on
c                 --- set zi
                  call twh_ntype(GLB_GET,imolty,iunit,itype)
                  ii = twh_getnbtype(itype,itype)
                  zi = twh_get_nbcoeff(10,ii)
c                 --- set zk
                  kunit = twh_get_ijben3(imolty,iunit,iben)
                  call twh_ntype(GLB_GET,imolty,kunit,ktype)
                  kk = twh_getnbtype(ktype,ktype)
                  zk = twh_get_nbcoeff(10,kk)
c                 --- set theta naught
                  itemp = twh_get_itben(imolty,iunit,iben)
                  theta = twh_get_bencoeff(itemp,0)
c                 --- get the equilibrium bond lengths
                  itemp = twh_get_ibtben1(imolty,iunit,iben)
                  rij = twh_get_vibcoeff(itemp,0)
                  itemp = twh_get_ibtben2(imolty,iunit,iben)
                  rjk = twh_get_vibcoeff(itemp,0)
c                 --- determine the rik squared value
                  rik = dsqrt(rij*rij + rjk*rjk 
     &                 - 2.0d0*rij*rjk*dcos(theta))
c                 --- set benprefact
                  dvalue = CON_KCALTOK*(664.12d0)
     &                 *(zi*zk/(rik**5))*( 
     &                 3.0d0*rij*rjk*(1.0d0 - dcos(theta)**2)
     &                 - rik*rik*dcos(theta))
                  if ( twh_get_anglestyle(atype)
     &                 .eq. ANG_STYLE_FOURIER ) then
c                    --- divide out by n^2
                     itemp = twh_get_itben(imolty,iunit,iben)
                     dvalue = dvalue/(twh_get_bencoeff(itemp,1)**2)
                  endif
                  call twh_benprefact(GLB_SET,imolty,iunit,iben,dvalue)
               endif
            enddo
         enddo

         do iunit = 1,nunit
c           --- check torsions
            do itor = 1, twh_get_intor(imolty,iunit)
               itest = twh_get_ijtor2(imolty,iunit,itor)
               itestb = twh_get_ijtor3(imolty,iunit,itor)
               itestc = twh_get_ijtor4(imolty,iunit,itor)
c              --- check that atom #2 is bonded to this atom
               lfound = .false.
               do ivib = 1,twh_get_invib(imolty,iunit)
                  if ( twh_get_ijvib(imolty,iunit,ivib) .eq. itest 
     &                 ) then
                     lfound = .true.
                     itemp = twh_get_itvib(imolty,iunit,ivib) 
                     call twh_ibttor1(GLB_SET,imolty,iunit,itor,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 1-2 connected in torsion are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #3 is bonded to atom #2
               lfound = .false.
               do ivib = 1,twh_get_invib(imolty,itest)
                  if ( twh_get_ijvib(imolty,itest,ivib) .eq. itestb 
     &                 ) then
                     lfound = .true.
                     itemp = twh_get_itvib(imolty,itest,ivib)
                     call twh_ibttor2(GLB_SET,imolty,iunit,itor,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 2-3 connected in torsion are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #4 is bonded to atom #3
               lfound = .false.
               do ivib = 1,twh_get_invib(imolty,itestb)
                  if ( twh_get_ijvib(imolty,itestb,ivib) .eq. itestc 
     &                 ) then
                     lfound = .true.
                     itemp = twh_get_itvib(imolty,itestb,ivib) 
                     call twh_ibttor3(GLB_SET,imolty,iunit,itor,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 3-4 connected in torsion are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that the 1-2-3 bend exists
               lfound = .false.
               do iben = 1, twh_get_inben(imolty,iunit)
                  if ( twh_get_ijben2(imolty,iunit,iben) .eq. itest 
     &                 .and.
     &                 twh_get_ijben3(imolty,iunit,iben) .eq. itestb 
     &                 ) then
                     lfound = .true.
                     itemp = twh_get_itben(imolty,iunit,iben)
                     call twh_iattor1(GLB_SET,imolty,iunit,itor,itemp) 
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 1-2-3 connected in torsion are not '
     &                 ,'connected by a bend'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that the 2-3-4 bend exists
               lfound = .false.
               do iben = 1, twh_get_inben(imolty,itest)
                  if ( twh_get_ijben2(imolty,itest,iben) .eq. itestb 
     &                 .and.
     &                 twh_get_ijben3(imolty,itest,iben) .eq. itestc 
     &                 ) then
                     lfound = .true.
                     itemp = twh_get_itben(imolty,itest,iben)
                     call twh_iattor2(GLB_SET,imolty,iunit,itor,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 2-3-4 connected in torsion are not '
     &                 ,'connected by a bend'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif

c              --- check that each other atom also contains this 
c              --- interaction
               lfound = .false.
               do jtor = 1, twh_get_intor(imolty,itestc)
                  if ( twh_get_ijtor4(imolty,itestc,jtor) .eq. iunit 
     &                 .and. 
     &                 twh_get_ijtor3(imolty,itestc,jtor) .eq. itest
     &                 .and.
     &                 twh_get_ijtor2(imolty,itestc,jtor) .eq. itestb 
     &                 ) then
                     if ( abs(twh_get_ittor(imolty,itestc,jtor)) .ne. 
     &                    abs(twh_get_ittor(imolty,iunit,itor)) ) then
                        write(6,*) 'Problem: imolty,iunit,itor',imolty
     &                       ,iunit,itor
                        write(6,*) 'forward atoms ',iunit,itest,itestb
     &                       ,itestc
                        write(6,*) 'reverse atoms ',itestc
     &                       ,twh_get_ijtor2(imolty,itestc,jtor)
     &                       ,twh_get_ijtor3(imolty,itestc,jtor)
     &                       ,twh_get_ijtor4(imolty,itestc,jtor)
                        write(6,*) 'iunit itor forward type '
     &                       ,twh_get_ittor(imolty,iunit,itor)
                        write(6,*) 'itestb jtor reverse type '
     &                       ,twh_get_ittor(imolty,itestc,jtor)
                        if ( lintraterm ) then
                           lfailure = .true.
                           return
                        endif
                     endif
                     lfound = .true.
                  endif
               enddo
               if ( .not. lfound ) then 
                  write(6,*) 'Problem: imolty,iunit,itor',imolty
     &                 ,iunit,itor
                  write(6,*) 'There is no symmetric torsion'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
         enddo

c        --- setup the angletypes for the angle-angle terms
         do iunit = 1,nunit
            do iaa = 1,twh_get_inaa(imolty,iunit)
               itest = twh_get_ijaa0(imolty,iunit,iaa)
               itestb = twh_get_ijaa1(imolty,iunit,iaa)
               itestc = twh_get_ijaa2(imolty,iunit,iaa)
               do iben = 1,twh_get_inben(imolty,itest)
                  if ( iunit .eq. twh_get_ijben2(imolty,itest,iben) 
     &                 ) then
                     if (itestb .eq. twh_get_ijben3(imolty,itest,iben)
     &                    ) then
                        itemp = twh_get_itben(imolty,itest,iben)
                        call twh_iataa1(GLB_SET,imolty,iunit,iaa,itemp)
                     elseif( itestc .eq. 
     &                       twh_get_ijben3(imolty,itest,iben)) then
                        itemp = twh_get_itben(imolty,itest,iben)
                        call twh_iataa2(GLB_SET,imolty,iunit,iaa,itemp)
                     endif
                  endif
               enddo
            enddo
         enddo

c        --- check improper torsions
         do iunit = 1,nunit
            do iimp = 1, twh_get_inimprop(imolty,iunit)
c              --- check that atom #2 is bonded to first atom
               itest = twh_get_ijimprop2(imolty,iunit,iimp)
               lfound = .false.
               do ivib = 1,twh_get_invib(imolty,iunit)
                  if ( twh_get_ijvib(imolty,iunit,ivib) .eq. itest ) 
     &                 lfound = .true.
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijimprop2,iimp'
     &                 ,imolty,iunit,itest,iimp
                  write(6,*) 'Atoms 1-2 connected in improper are '
     &                 ,'not connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #3 is bonded to first atom
               itestb = twh_get_ijimprop3(imolty,iunit,iimp)
               lfound = .false.
               do ivib = 1,twh_get_invib(imolty,iunit)
                  if ( twh_get_ijvib(imolty,iunit,ivib) .eq. itestb ) 
     &                 lfound = .true.
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijimprop3,iimp'
     &                 ,imolty,iunit,itestb,iimp
                  write(6,*) 'Atoms 1-3 connected in improper are '
     &                 ,'not connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #4 is bonded to first atom
               itestc = twh_get_ijimprop4(imolty,iunit,iimp)
               lfound = .false.
               do ivib = 1,twh_get_invib(imolty,iunit)
                  if ( twh_get_ijvib(imolty,iunit,ivib) .eq. itestc ) 
     &                 lfound = .true.
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijimprop4,iimp'
     &                 ,imolty,iunit,itestc,iimp
                  write(6,*) 'Atoms 1-4 connected in improper are '
     &                 ,'not connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
         enddo
c        --- check special one-five
         do iunit = 1,nunit
            do iof = 1,twh_get_inof(imolty,iunit)
               matom = twh_get_ijof5(imolty,iunit,iof)
               lfound = .false.
               do jof = 1,twh_get_inof(imolty,matom)
                  if ( twh_get_ijof5(imolty,matom,jof) .eq. iunit ) 
     &                 lfound = .true.
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'CHECKSTRUC: non-symetric one-five'
                  write(6,*) 'imolty,iunit,iof',imolty,iunit,iof
                  write(6,*) 'has no return interaction from ijof5'
     &                 ,twh_get_ijof5(imolty,iunit,iof)
                  lfailure = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- need some special logic for torsions that have 1-4 nonbonded
c     --- interactions.  In certain cases (6-membered rings) an atom
c     --- can have multiple torsions which end at the same 4th bead
c     --- In this case we need to scale the nonbonded interactions to 
c     --- avoid counting the nonbond terms multiple times
c     --- for other ring structures it is possible that the 1-4 is also
c     --- a 1-3 or 1-2 (smaller rings) and we want to disable the
c     --- nonbonded calculation
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
c           --- initialize the array of cross terms
            do junit = 1,nunit
               torcross(junit) = 0
            enddo
            call twh_intor(GLB_GET,imolty,iunit,numtor)
            do itor = 1,numtor
c              --- check to see if there are nonbonded terms for this
c              --- torsion
               call twh_ittor(GLB_GET,imolty,iunit,itor,itype)
               call twh_loftor(GLB_GET,iabs(itype),lonefourtor)
               if ( lonefourtor ) then
                  torcross(twh_get_ijtor4(imolty,iunit,itor)) = 
     &                 torcross(twh_get_ijtor4(imolty,iunit,itor)) + 1
               endif
            enddo
c           --- now make sure we don't have bonded or angle interaction
c           --- with that 4th bead (possible in rings)
            do itor = 1,twh_get_intor(imolty,iunit)
               junit = twh_get_ijtor4(imolty,iunit,itor)
               lfound = .false.
               do ivib = 1,twh_get_invib(imolty,iunit)
                  if ( junit .eq. twh_get_ijvib(imolty,iunit,ivib) )
     &                 lfound = .true.
               enddo
               do iben = 1,twh_get_inben(imolty,iunit)
                  if ( junit .eq. twh_get_ijben3(imolty,iunit,iben) ) 
     &                 lfound = .true.
               enddo
c              --- store the scaling factor for each torsion
               if ( lfound .or. torcross(junit) .eq. 0 ) then
c                 --- we do not want to count these nonbonded interact
                  call twh_torofcode(GLB_SET,imolty,iunit,itor
     &                 ,CODE_ZERO)
               else
c                 --- use the torcross to compute scale factor
                  if ( torcross(junit) .gt. 1 ) then
                     dvalue = 1.0d0 / dble(torcross(junit))
                     call twh_torofcode(GLB_SET,imolty,iunit,itor
     &                    ,CODE_RESCALE)
                     call twh_itscale(GLB_SET,imolty,iunit,itor,dvalue)
                  endif
               endif
            enddo
         enddo
      enddo

c     --- check that the hybridization of the atoms is consistent with
c     --- the number of bonds.  Note that since united-atom force fields
c     --- are possible there may not be an exact match, this is more of
c     --- a less than or equal to senario
c     --- we mostly just print out warnings, but if we have too many 
c     --- bonds to an atom for its hybridization then we are in trouble
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            ivib = twh_get_invib(imolty,iunit)
            call twh_ntype(GLB_GET,imolty,iunit,itype)
            if ( twh_get_bondpatt(itype) .eq. 'sp3' ) then
               if ( ivib .gt. 4 ) then
                  write(6,*) 'too many bonds to an sp3 atom'
                  write(6,*) 'imolty,iunit,ivib',imolty,iunit,ivib
                  write(6,*) 'itype,bondpatt',itype
     &                 ,twh_get_bondpatt(itype)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            elseif ( twh_get_bondpatt(itype) .eq. 'sp2' ) then
               if ( ivib .gt. 3 ) then
                  write(6,*) 'too many bonds to an sp2 atom'
                  write(6,*) 'imolty,iunit,ivib',imolty,iunit,ivib
                  write(6,*) 'itype,bondpatt',itype
     &                 ,twh_get_bondpatt(itype)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            elseif ( twh_get_bondpatt(itype) .eq. 'sp' ) then
               if ( ivib .gt. 2 ) then
                  write(6,*) 'too many bonds to an sp atom'
                  write(6,*) 'imolty,iunit,ivib',imolty,iunit,ivib
                  write(6,*) 'itype,bondpatt',itype
     &                 ,twh_get_bondpatt(itype)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            elseif ( twh_get_bondpatt(itype) .eq. 's' ) then
               if ( ivib .gt. 1 ) then
                  write(6,*) 'too many bonds to an s atom'
                  write(6,*) 'imolty,iunit,ivib',imolty,iunit,ivib
                  write(6,*) 'itype,bondpatt',itype
     &                 ,twh_get_bondpatt(itype)
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            endif
         enddo
      enddo

c     --- check that the molecule is fully connected and there is no
c     --- subset of the molecule that is not connected to the rest.
      do imolty = 1, nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         call twh_checkconnect(imolty,nunit)
         if ( lfailure ) return
      enddo

      return
      end

      subroutine twh_checkconnect(imolty,nunit)
c     ******************************************************************
c     * checks the bonding pattern of a molecule to make sure it is    *
c     * fully connected.  Stores this information in lbondgraph        *
c     *                                                                *
c     * originally written 05-30-2013 by M.G. Martin                   *
c     * last modified 07-15-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer imolty,nunit
c     --- local variables
      logical lexist,lbondgraph
      integer workcount,ivib,numvib,iother,icount,icurrent,iunit
      integer worklist
      dimension worklist(nunit*NNBOND)

c     --- initialise lexist to all false
      call twh_lexist(GLB_INIT,nunit,.false.)
c     --- start with atom 1 and progress through the bonded stucture
      workcount = 1
      worklist(workcount) = 1

      do while ( workcount .gt. 0 ) 
         icurrent = worklist(1)
         call twh_lexist(GLB_GET,icurrent,lexist)
         if ( lexist ) then
c           --- this unit has already been found to exist and tested
c           --- nothing else to do here
         else
c           --- add all of the atoms bonded to icurrent to the worklist
            call twh_invib(GLB_GET,imolty,icurrent,numvib)
            do ivib = 1,numvib
               call twh_ijvib(GLB_GET,imolty,icurrent,ivib,iother)
               workcount = workcount + 1
               worklist(workcount) = iother
            enddo
c           --- set this atom to exist
            call twh_lexist(GLB_SET,icurrent,.true.)
         endif
c        --- decrement and shift the worklist
         do icount = 1,workcount - 1
            worklist(icount) = worklist(icount+1)
         enddo
         workcount = workcount - 1
      enddo

      lbondgraph = .true.
      do iunit = 1,nunit
         call twh_lexist(GLB_GET,iunit,lexist)
         if ( .not. lexist ) then
c           --- this atom was never found in the pattern
c           --- warn on first occurance
            if ( lbondgraph ) then
               write(6,'(a,1x,i4,1x,a,1x,i6,1x,a)')
     &              'Molecule type',imolty,'unit'
     &              ,iunit,' not connected via bonding graph to atom 1'
               write(6,'(2x,a,1x,a)')
     &              'therefore Monte Carlo moves will not be allowed'
     &              ,'to alter atom positions for this molecule type'
               lbondgraph = .false.
            endif
         endif
      enddo
c     --- set lbondgraph
      call twh_lbondgraph(GLB_SET,imolty,lbondgraph)

      return
      end

