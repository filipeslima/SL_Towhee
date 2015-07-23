#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2012 Marcus G. Martin                       *
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
      subroutine twh_mapconf(lfailure)
c     ******************************************************************
c     * this subroutine creates or reads a towhee_map file             *
c     *                                                                *
c     * originally created 02-2006 by O. Yazaydin                      *
c     * last modified 09-18-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#include "functions.h"
c     --- debug flag
#define DEBUG_MAPCONF 0
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- logical scalars
      logical loverlap,lcreatemap
c     --- logical arrays
      logical lnbterm
      dimension lnbterm(1)
c     --- integer scalars
      integer cubea,cubeb,cubec,testcord,globalcord,cubeindex
      integer engastyle,testchain,nchoi,iufrom,ntogrow,ibox,nmolty
      integer iunit,nchain,nunit,mapmolty,iaa,ibb,icc
      integer cubetotal
c     --- double precision scalars
      double precision arba,arbb,arbc,xcen,ycen,zcen
      double precision beta,cweight,rdist,cweighttotal
      double precision totalweight,dcube,sizea,sizeb,sizec
c     --- double precision arrays
      double precision rxp,ryp,rzp
      dimension rxp(1,1),ryp(1,1),rzp(1,1)
      double precision vctry,vtrext,vtrintra,vtrinter,vtrelect,vtrewald
      dimension  vctry(1),vtrext(1),vtrintra(1),vtrinter(1),vtrelect(1)
     &          ,vtrewald(1)
      double precision vqtry
      dimension vqtry(1)

c     --- retrieve constants
      call twh_nmolty(GLB_GET,nmolty)
      call twh_nchain(GLB_GET,nchain)
      call twh_lcreatemap(GLB_GET,lcreatemap)

      if (lcreatemap) then
c        --- retrieve the map molecule type
         call twh_mapmolty(GLB_GET,mapmolty)
c        --- get the number of atoms in this molecule type
         call twh_nunit(GLB_GET,mapmolty,nunit)
c        --- set testflag to use the single molecule coordinates
         testcord = CRD_SINGLE
         globalcord = CRD_REAL
         engastyle = 1
         testchain = nchain + 1
         call twh_moltyp(GLB_SET,testchain,mapmolty)
         nchoi = 1
         iufrom = 0
         ntogrow = 1
         rdist = 0.0d0
c        --- set lexist to false for all atoms in this molecule type
         call twh_lexist(GLB_INIT,nunit,.false.)
c        ---set some necessary variables to create the map
         ibox = 1
c        ---initialize weight variables
         call twh_cubelet(GLB_GET,cubea,cubeb,cubec)
c        --- compute subcube sizes
         sizea = 1.0d0/dble(cubea)
         sizeb = 1.0d0/dble(cubeb)
         sizec = 1.0d0/dble(cubec)
c        --- allocate memory for total cubelets
         cubetotal = cubea*cubeb*cubec
         cweighttotal = 0.0d0
         call twh_cubeletweight(GLB_ALLOC,cubetotal,cweighttotal)

         write(6,'(a40)') 'MAPCONF: creating the energy biasing map'
c        --- compute the average Boltzmann weight for each cubelet
         totalweight = 0.0d0
         do iaa = 1, cubea
            do ibb = 1, cubeb
               do icc = 1, cubec
c                 --- there are many ways that could be used to generate
c                 --- the biasing map.  For now, we are taking the 
c                 --- first bead insertion Boltzmann weight averaged
c                 --- over all possible first atoms in the molecule
                  arba = (dble(iaa-1)+0.5d0)*sizea
                  arbb = (dble(ibb-1)+0.5d0)*sizeb
                  arbc = (dble(icc-1)+0.5d0)*sizec
                  call twh_arbtocart(ibox,arba,arbb,arbc
     &                 ,xcen,ycen,zcen)
c                 --- store this in the r*p arrays
                  rxp(1,1) = xcen
                  ryp(1,1) = ycen
                  rzp(1,1) = zcen
                  cweighttotal = 0.0d0
                  do iunit = 1,nunit
c                    --- set the active atom to the current unit
                     call twh_glist(GLB_SET,1,iunit)
                     lnbterm(1) = .false.
c                    --- compute the insertion energy into the box
                     call twh_engatom(lfailure,testcord,globalcord
     &                    ,engastyle,loverlap,testchain,mapmolty
     &                    ,ibox,nchoi,iufrom,ntogrow,rdist,lnbterm
     &                    ,1,1,rxp,ryp,rzp,vctry,vtrext,vtrintra
     &                    ,vtrinter,vtrelect,vtrewald,vqtry)
                     if ( lfailure ) return
c                    --- compute the energy for each cube
                     if ( loverlap ) then
c                       --- hard overlap, set to zero weight
                        cweight = 0.0d0
                     else
c                       --- compute weight
                        call twh_beta(GLB_GET,beta)
                        cweight = twh_expon(-beta*vctry(1)) 
                     endif
                     cweighttotal = cweighttotal + cweight
                  enddo
c                 --- store the unnormalized weight for this subcube
                  call twh_cubelet_index(iaa,ibb,icc,cubeindex)
                  call twh_cubeletweight(GLB_SET,cubeindex,cweighttotal)
c                 --- increment the totalweight
                  totalweight = totalweight + cweighttotal
               enddo
            enddo
         enddo 

         if ( totalweight .lt. SMALLEST ) then
            write(6,*) 'MAPCONF: total Rosenbluth weight low:'
     &           ,totalweight
            write(6,*) 'using a uniform distribution of cubes instead'
            cweight = 1.0d0
            totalweight = 0.0d0
            do iaa = 1,cubea
               do ibb = 1,cubeb
                  do icc = 1,cubec
                     call twh_cubelet_index(iaa,ibb,icc,cubeindex)
                     call twh_cubeletweight(GLB_SET,cubeindex,cweight)
                     totalweight = totalweight + cweight
                  enddo
               enddo
            enddo
         endif

c        --- normallize the distribution
         dcube = dble(cubea)*dble(cubeb)*dble(cubec)
         cweighttotal = 0.0d0
         do iaa = 1, cubea
            do ibb = 1, cubeb
               do icc = 1, cubec
                  call twh_cubelet_index(iaa,ibb,icc,cubeindex)
                  call twh_cubeletweight(GLB_GET,cubeindex,cweight)
                  cweight = cweight/totalweight
c                 --- mix a 99% chance of energybias with a 1% chance
c                 --- of randomly selecting a cubelet to make sure the
c                 --- bias is positive for all cubelets
                  cweight = 0.99d0*cweight + 0.01d0/dcube
                  call twh_cubeletweight(GLB_SET,cubeindex,cweight)
c                 --- add up weights to double check
                  cweighttotal = cweighttotal + cweight
               enddo
            enddo
         enddo

         write(6,'(a22,1x,f12.10)') 'MAPCONF: cweighttotal:'
     &        ,cweighttotal
         if ( dabs(cweighttotal-1.0d0) .gt. 1.0d-5 ) then
            write(6,*) 'Normalization failed'
            lfailure = .true.
            return
         endif
c        --- write the map information to the file 
         call twh_openfile(lfailure,68,'towhee_map','DEFAULT',
     &        'UNKNOWN','DEFAULT')
         if ( lfailure ) return

         write(68,*) cubea, cubeb, cubec
         do iaa = 1, cubea
            do ibb = 1, cubeb
               do icc = 1, cubec
                  call twh_cubelet_index(iaa,ibb,icc,cubeindex)
                  call twh_cubeletweight(GLB_GET,cubeindex,cweight)
                  write(68,*) cweight
               enddo
            enddo
         enddo
         rewind(68)
         close(68)

         write(6,'(a36)') 'MAPCONF: new towhee_map file created'
      else
c        --- read map info from file
         call twh_openfile(lfailure,68,'towhee_map','DEFAULT','OLD'
     &        ,'DEFAULT')
         if ( lfailure ) return
         read(68,*) cubea, cubeb, cubec
         call twh_cubelet(GLB_SET,cubea,cubeb,cubec)
c        --- allocate memory for total cubelets
         cubetotal = cubea*cubeb*cubec
         cweighttotal = 0.0d0
         call twh_cubeletweight(GLB_ALLOC,cubetotal,cweighttotal)

         do iaa = 1, cubea
            do ibb = 1, cubeb
               do icc = 1, cubec
                  read(68,*) cweight
                  call twh_cubelet_index(iaa,ibb,icc,cubeindex)
                  call twh_cubeletweight(GLB_SET,cubeindex,cweight)
               enddo
            enddo
         enddo

         rewind(68)
         close(68)
         write (6,'(a34)') 'MAPCONF: towhee_map read from file'
      endif
      return
      end

      subroutine twh_cubelet_index(ixx,iyy,izz,index)
c     ******************************************************************
c     * computes the cublet index from the three components            *
c     *                                                                *
c     * originally written 12-05-2008 by M.G. Martin                   *
c     * last modified 12-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ixx,iyy,izz,index
c     --- local variables
      integer xsize,ysize,zsize

c     --- convert from 3-D to single dimension array
      call twh_cubelet(GLB_GET,xsize,ysize,zsize)
      index = (ixx-1)*(ysize*zsize) + (iyy-1)*(zsize) + izz

      return
      end
