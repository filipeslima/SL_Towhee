#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2014 Marcus G. Martin                       *
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
      subroutine twh_ffreadcharmmfile(lfailure,nbtype,botype,antype
     &     ,ditype,imtype,ffname)
c     ******************************************************************
c     * reads in a charmm-type force field file and reads it into      *
c     * corresponding Towhee data structures so that it can create     *
c     * a Towhee force field file                                      *
c     *                                                                *
c     * written by M.G. Martin and C.D. Lorenz 02-18-03                *
c     * last modifed 07-21-2014 by M.G. Martin                         *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nbtype,botype,antype,ditype,imtype
      character*(FFNAMELEN) ffname
c     --- local variables
      dimension dataline(3,7)
      integer dataline,i,iline,idone,j,jline,jdone
      integer numlines,iflag,n,ii
      integer nvslocal,nbslocal,ntslocal
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
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)

      double precision fdummy,epsilon,rmin2
      double precision kb,bo,ktheta,thetao,kd,phio,psio
      double precision ki,onepi
      double precision masslocal
      character*2 elelocal 

      character*3 catom1,catom2,catom3,catom4
      character*10 charscan 
      character*30 filename
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
      integer nislocal 
      character*(FFNAMELEN) ninlocal 
      dimension ninlocal(MAXSAMEIMP,MAXIMPNAME)

      iflag = 0
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

      call twh_constant_pi(GLB_GET,onepi)

      write(6,*) 'Enter the filename of the charmm file to translate'
      read(5,*) filename
      write(6,*) 'Enter a name for your forcefield'
      read(5,*) ffname

      do i = 1,3
         do j = 1,7
            dataline(i,j) = 0
         enddo
      enddo

      numlines = 0 
      iline = 0
      idone = 1
      call twh_openfile(lfailure,34,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
      do i = 1,10000
        if (idone .gt. 0) then
           read(34,10) charscan
           numlines = numlines + 1           
           if ((charscan(1:1) .ne. "*") .and. (charscan(1:1) .ne. "!") 
     &         .and. (charscan(1:2) .ne. "  ")) then
               if (charscan .eq. "BONDS") then
                  iline = iline + 1
                  iflag = 2
               else if (charscan .eq. "ANGLES") then
                  iline = iline + 1
                  iflag = 3
               else if (charscan .eq. "DIHEDRAL") then
                  iline = iline + 1
                  iflag = 4
               else if (charscan .eq. "IMPROPER") then
                  iline = iline + 1
                  iflag = 5 
               else if (charscan .eq. "NONBONDED") then
                  iline = iline + 1
                  iflag = 1
               else if (charscan .eq. "NBFIX") then
                  iline = iline + 1
                  iflag = 6
               else if (charscan .eq. "END") then
                  idone = -1
               else 
                  if (dataline(1,iline) .eq. 0) then
                     dataline(1,iline) = iflag
                  endif
                  if (dataline(2,iline) .eq. 0) then
                     dataline(2,iline) = numlines
                  endif
                  if (i .gt. dataline(3,iline)) then
                    dataline(3,iline) = numlines
                  endif
               endif
           endif
        endif
      enddo

     
  10  format(a10)
      rewind 34

      jdone = 1
      jline = 1

      ii = 0

      do j = 1,numlines
         if (j .lt. dataline(2,jline)) then
            read(34,*)
         else 
            if (dataline(1,jline) .eq. 1) then
c              *** read in data for nonbond forcefield terms
               nbtype = nbtype + 1
               read(34,*) catom1,fdummy,epsilon,rmin2
               ii = twh_getnbtype(nbtype,nbtype)
               nbclocal(1) = rmin2 * (2.0d0)/(2.0d0**(1.0/6.0))
               nbclocal(2) = epsilon * (CON_KCALTOK)
               nbclocal(3) = nbclocal(1)
               nbclocal(4) = nbclocal(2)
               call twh_set_nbcoeff(ii,nbclocal)
               call twh_nbff(GLB_SET,nbtype,ffname)
               nbnlocal = catom1
               cbnlocal = nbnlocal
               canlocal = nbnlocal
               ctnlocal = nbnlocal
               call twh_set_names(nbtype,nbnlocal,cbnlocal,canlocal
     &              ,ctnlocal)
               if (catom1(1:1) .eq. "C") then
                  masslocal = 12.011d0
                  elelocal = ' C'
               else if (catom1(1:1) .eq. "N") then
                  masslocal = 14.007d0
                  elelocal = ' N'
               else if (catom1(1:1) .eq. "O") then
                  masslocal = 15.999d0
                  elelocal = ' O'
               else if (catom1(1:1) .eq. "S") then
                  masslocal = 32.066d0
                  elelocal = ' S'
               else if (catom1(1:1) .eq. "H") then
                  masslocal = 1.008d0
                  elelocal = ' H'
               else
                  write(6,*) "ERROR: Unknown element type",
     &              catom1(1:1)
                  lfailure = .true.
                  return
               endif
               call twh_elementname(GLB_SET,nbtype,elelocal) 
               call twh_mass(GLB_SET,nbtype,masslocal) 
            else if (dataline(1,jline) .eq. 2) then
c              *** read in data for bond forcefield terms
               read(34,*) catom1,catom2,kb,bo
               botype = botype + 1
               call twh_vibstyle(GLB_SET,botype,BOND_STYLE_HARMONIC)
               vclocal(0) = bo
               vclocal(1) = kb*(CON_KCALTOK)
               call twh_set_vibcoeffs(botype,vclocal)
               call twh_vibff(GLB_SET,botype,ffname)
               nvslocal = 1
               vnlocal(1,1) = catom1
               vnlocal(1,2) = catom2
               call twh_set_vibnames(botype,nvslocal,vnlocal)
            else if (dataline(1,jline) .eq. 3) then
c              *** read in data for angle forcefield terms
               read(34,*) catom1,catom2,catom3,ktheta,thetao
               antype = antype + 1
               call twh_anglestyle(GLB_SET,antype,ANG_STYLE_HARMONIC)
               bclocal(0) = thetao
               bclocal(1) = ktheta * (CON_KCALTOK)
               call twh_set_bencoeff(antype,bclocal)
               call twh_bendff(GLB_SET,antype,ffname)
               nbslocal = 1
               bnlocal(1,1) = catom1
               bnlocal(1,2) = catom2
               bnlocal(1,3) = catom3
               call twh_set_bendnames(antype,nbslocal,bnlocal)
            else if (dataline(1,jline) .eq. 4) then
c           *** read in data for dihedral forcefield terms
                read(34,*) catom1,catom2,catom3,catom4,kd,n,phio
                ditype = ditype + 1
                call twh_torstyle(GLB_SET,ditype,3)
                call twh_loftor(GLB_SET,ditype,.true.)
                tclocal(-1) = 1.0d0
                call twh_ntorloop(GLB_SET,ditype,1)
                tclocal(1) = kd*(CON_KCALTOK)
                tclocal(2) = n * 1.0d0
                tclocal(3) = phio * (onepi/180.0d0)
                call twh_set_torcoeff(ditype,tclocal)
                call twh_torff(GLB_SET,ditype,ffname)
                ntslocal = 1
                tnlocal(1,1) = catom1
                tnlocal(1,2) = catom2
                tnlocal(1,3) = catom3
                tnlocal(1,4) = catom4
                call twh_set_tornames(ditype,ntslocal,tnlocal)
            else if (dataline(1,jline) .eq. 5) then
c           *** read in data for improper forcefield terms
                read(34,*) catom1,catom2,catom3,catom4,ki,n,psio
                imtype = imtype + 1
                call twh_impform(GLB_SET,imtype,2)
                call twh_impstyle(GLB_SET,imtype,1)
                iclocal(0) = ki * (CON_KCALTOK)
                iclocal(1) = psio * (onepi/180.0d0)
                call twh_set_impcoeff(imtype,iclocal)
                call twh_impff(GLB_SET,imtype,ffname)
                nislocal = 1
                ninlocal(1,1) = catom1
                ninlocal(1,2) = catom2
                ninlocal(1,3) = catom3
                ninlocal(1,4) = catom4
                call twh_set_impnames(imtype,nislocal,ninlocal)
            endif
            if (j .eq. dataline(3,jline)) then
               jline = jline + 1
            endif   
         endif
      enddo                        
      close(34)

      return
      end

