#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2011 Marcus G. Martin                       *
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
      subroutine twh_initconf(lfailure)
c     ******************************************************************
c     * this subroutine sets up a simple initial conformation for      *
c     * the system.                                                    *
c     * initstyle(ibox,imolty) must be one of the following            *
c     * 'full cbmc': create template via cbmc from scratch             *
c     * 'template': read molecule template from towhee_template        *
c     * 'coords': read all coords from towhee_coords                   *
c     * 'nanotube': towhee_nanotube molecule template                  *
c     * 'helix cbmc': create towhee_template by placing some atoms     *
c     *               using buildhelix and then grow the rest with cbmc*
c     * 'partial cbmc': read in atom coordinates from towhee_partial   *
c     *                 and fill in any missing atoms with CBMC        *
c     * 'unit cell': read in a template of molecule atomic positions   *
c     *              from towhee_cell and duplicate the cell to create *
c     *              initial file                                      *
c     * if initstyle is negative then the center of mass of the        *
c     *    molecule is placed in the center of the box                 *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modified 08-16-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- set the file identities
#define IFILE_TEMPLATE 78
#define IFILE_COORDS 79
#define IFILE_NANO 80
#define IFILE_PARTIAL 81
#define IFILE_UNITCELL 82
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- character scalars
      character*15 idname
      character*20 initboxtype,initstyle,teststyle,initlattice
c     --- logical scalars
      logical lopentemplate,lopencoords,lopenpartial,lallow
      logical lterm,lgrowcbmc,lfound,lcoulomb,lewald
c     --- integer scalars
      integer idim,jdim,istart,iend,iadjust,itest,nchain,kmax
      integer ic,jc,kc,intemp,idum,nchbox,jbox,ensemble
      integer ibox,unitc,imolty,ichain,isolvtemp
      integer iunit,isteps,isafety,iutry,movetype,nmtype
      integer fieldtemp,ntest,testcord,idealchain,uc_count,nmolty,nunit
      integer inix,iniy,iniz,initmol,globalpos,itemp,numboxes
      integer atomcount,boxtest
c     --- integer arrays
      integer mcmt,mcmtma
      dimension mcmt(NTMAX,MAXBOX),mcmtma(NTMAX,MAXBOX)
      integer cell_initmol
      dimension cell_initmol(NTMAX)
c     --- double precision scalars
      double precision arba,arbb,arbc,xtemp,ytemp,ztemp,dpscale
      double precision xcen,ycen,zcen,xoffset,yoffset,zoffset
      double precision xcmi,ycmi,zcmi,calp,boxvolume
c     --- double precision arrays
      double precision ux,uy,uz
      dimension ux(MAXBOX),uy(MAXBOX),uz(MAXBOX)
      double precision samx,samy,samz
      dimension samx(NTMAX,NUMAX),samy(NTMAX,NUMAX),samz(NTMAX,NUMAX)
      double precision cell_hmatrix
      dimension cell_hmatrix(3,3)
      double precision dimtemp
      dimension dimtemp(3)
      double precision vrostot,vrosvib,vrosang,vrostor,vrosext
      dimension vrostot(1),vrosvib(1),vrosang(1),vrostor(1),vrosext(1)
      double precision vrosintra,vrosinter,vroselect,vrosewald
      dimension vrosintra(1),vrosinter(1),vroselect(1),vrosewald(1)
      double precision wrosen
      dimension wrosen(2)

c      write(6,*) 'starting INITCONF'

c     --- retrieve constants
      call twh_nmolty(GLB_GET,nmolty)
      call twh_numboxes(GLB_GET,numboxes)

c     --- set idname
      idname = 'initconf'
c     --- set testflag to use the single molecule coordinates
      testcord = CRD_SINGLE
      call twh_random_allow_restart(GLB_GET,lallow)
      if ( lallow ) then
c        --- initialize the random number generator on a singe integer
         call twh_init_random_single(lfailure)
         if ( lfailure ) return
      endif

      call twh_initboxtype(GLB_GET,initboxtype)
      if ( initboxtype .eq. 'unit cell' ) then
         call twh_ini_read_unitcell(lfailure,cell_hmatrix, cell_initmol)
         if ( lfailure ) return
      endif

c     --- invert and safetycheck the hmatrix
      do ibox = 1,numboxes
c        --- check this hmatrix
         call twh_checkhmatrix(ibox,.true.,lfailure)
         if ( lfailure ) then
            write(6,*) 'INITCONF: problem with unit cell hmatrix'
            return
         endif

c        --- invert this matrix
         call twh_inverthmatrix(ibox)
         do idim = 1,3
            do jdim = 1,3
               call twh_hinverse(GLB_GET,ibox,idim,jdim,dimtemp(jdim))
            enddo
            write(6,'(a18,i5,1x,i1,3(1x,f10.5))') 
     &           'Box idim hinverse:'
     &           ,ibox,idim,(dimtemp(jdim),jdim=1,3)
         enddo
c        --- make sure the volumes are positive
         call twh_boxvolume(GLB_GET,ibox,boxvolume)
         if ( boxvolume .lt. 0.0d0 ) then
            write(6,*) 'INITCONF: Box ',ibox
     &           ,' has a negative volume of ',boxvolume
            write(6,*) 'this means you have set up a left-handed '
     &           ,'coordinate system'
            lfailure = .true.
            return
         endif
      enddo

c     --- set nchbox
      do imolty = 1,nmolty
         do ibox = 1,numboxes
            call twh_initmol(GLB_GET,ibox,imolty,initmol)
            call twh_nchbox(GLB_INCR,ibox,initmol)
         enddo
      enddo

      if ( initboxtype .ne. 'unit cell' ) then
c        --- make sure the cell of the molecules is large enough to 
c        --- hold them
         do ibox = 1,numboxes
            call twh_inixyz(GLB_GET,ibox,inix,iniy,iniz)
            unitc = inix*iniy*iniz
            call twh_nchbox(GLB_GET,ibox,nchbox)
            if ( nchbox .gt. unitc ) then
               write(6,*) 'INITCONF: unit cell too small in box',ibox
               write(6,*) 'nchbox,unitc',nchbox,unitc
               lfailure = .true.
               return
            endif
         enddo
c        --- calculation of unit cell dimensions
         do ibox = 1,numboxes
            call twh_inixyz(GLB_GET,ibox,inix,iniy,iniz)
            ux(ibox) = 1.0d0 / dble(inix) 
            uy(ibox) = 1.0d0 / dble(iniy)
            uz(ibox) = 1.0d0 / dble(iniz)
         enddo
      endif

c     --- set mcmt to the specified initial values
      do imolty = 1,nmolty
         do ibox = 1,numboxes
            call twh_initmol(GLB_GET,ibox,imolty,initmol)
            mcmt(imolty,ibox) = initmol
         enddo
      enddo
c     --- set mcmtma to be the total number of molecules that type 
c     --- and below
      do ibox = 1,numboxes
         mcmtma(1,ibox) = mcmt(1,ibox)
         do imolty = 2, nmolty
            mcmtma(imolty,ibox) = mcmtma(imolty-1,ibox) 
     &           + mcmt(imolty,ibox)
         enddo
      enddo

      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_lewald(GLB_GET,lewald)
      if ( lcoulomb ) then
c        --- set ewald parameters: 
         do ibox = 1, numboxes
            call twh_resetcoulomb(lfailure,ibox)
            if ( lfailure ) return
            if ( lewald ) then
               call twh_ewald_calp(GLB_GET,ibox,calp)
               write(6,41) 'Box:',ibox,'Initial calp:',calp
 41            format(a4,1x,i5,1x,a13,1x,f10.5)
               call twh_ewald_kmax(GLB_GET,ibox,kmax)
               write(6,'(a4,1x,i5,1x,a13,1x,i5)') 
     &              'Box:',ibox,'Initial kmax:',kmax
            endif
         enddo
      endif

c     --- initialize nboxi
      call twh_nchain(GLB_GET,nchain)
      call twh_nboxi(GLB_INIT,nchain,0)

c     --- initialize atomcount
      atomcount = 0
      globalpos = 3*atomcount+1
      call twh_globalpos(GLB_SET,1,globalpos)
      
c     --- logicals to see if we have opened files
      lopentemplate = .false.
      lopencoords = .false.
      lopenpartial = .false.

c     --- set ichain to 0
      ichain = 0
      do 102 ibox = 1,numboxes
c        --- determine the type of template for each of the molecules 
c        --- and whether we can reuse this template from a previous 
c        --- box
         do imolty = 1, nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            call twh_initmol(GLB_GET,ibox,imolty,initmol)
            call twh_initstyle(GLB_GET,ibox,imolty,initstyle)
            if ( initmol .ne. 0 ) then
c              --- see if we already have this template from the 
c              --- previous box
               if ( ibox .gt. 1 ) then
                  jbox = ibox - 1
                  call twh_initmol(GLB_GET,jbox,imolty,itest)
                  call twh_initstyle(GLB_GET,jbox,imolty,teststyle)
                  if ( (itest .ne. 0) .and.
     &                 (teststyle .eq. initstyle) ) then
c                     --- no need to do this again, just use the 
c                     --- previous box template
                     goto 100
                  endif
               endif
               lgrowcbmc = .false.
               if ( initstyle .eq. 'full cbmc' ) then
                  lgrowcbmc = .true.
c                 --- grow the entire molecule using CBMC
c                 --- set the first atom to the origin and grow
c                 --- from there
                  xtemp = 0.0d0
                  ytemp = 0.0d0
                  ztemp = 0.0d0
                  call twh_coordinates(GLB_SET,testcord,ichain,1
     &                 ,xtemp,ytemp,ztemp)
                  movetype = 2
                  iutry = 1
               elseif ( initstyle .eq. 'template' ) then
c                 --- see if the towhee_template file needs to be opened
                  if ( .not. lopentemplate ) then
                     call twh_openfile(lfailure,IFILE_TEMPLATE
     &                    ,'towhee_template','DEFAULT','OLD','DEFAULT')
                     if ( lfailure ) return
                     lopentemplate = .true.
                  endif
c                 --- read sample structure from towhee_template
                  read(IFILE_TEMPLATE,*)
                  do iunit = 1, nunit
                     read(IFILE_TEMPLATE,*) samx(imolty,iunit)
     &                    ,samy(imolty,iunit),samz(imolty,iunit)
                  enddo
               elseif ( initstyle .eq. 'coords' ) then
c                 --- see if we need to open the coord file
                  if ( .not. lopencoords ) then
                     call twh_openfile(lfailure,IFILE_COORDS
     &                    ,'towhee_coords','DEFAULT','OLD','DEFAULT')
                     if ( lfailure ) return
                     lopencoords = .true.
                  endif
               elseif ( initstyle .eq. 'nanotube' ) then
c                 --- open the nanotube file
                  call twh_openfile(lfailure,IFILE_NANO
     &                 ,'towhee_nanotube','DEFAULT', 'OLD','DEFAULT')
                  if ( lfailure ) return
c                 --- read sample structure from towhee_nanotube
                  read(IFILE_NANO,*)
                  do iunit = 1, nunit
                     read(IFILE_NANO,*) samx(imolty,iunit), 
     &                     samy(imolty,iunit), samz(imolty,iunit)
                  enddo
c                 --- close the nanotube file
                  close(IFILE_NANO)
               elseif ( initstyle .eq. 'helix cbmc' ) then
                  lgrowcbmc = .true.
c                 --- create a helix using the keyatom
                  call twh_buildhelix(lfailure,imolty,testcord)
                  if ( lfailure ) return
c                 --- fill in the rest of the atoms using CBMC
                  iutry = 0
                  movetype = 5
               elseif ( initstyle .eq. 'partial cbmc') then
                  lgrowcbmc = .true.
                  if ( .not. lopenpartial ) then
                     call twh_openfile(lfailure,IFILE_PARTIAL
     &                    ,'towhee_partial','DEFAULT','OLD','DEFAULT')
                     if ( lfailure ) return
                     lopenpartial = .true.
                  endif
c                 --- seed with available atoms and then grow the rest
                  call twh_buildpartial(lfailure,imolty,testcord
     &                 ,IFILE_PARTIAL)
                  if ( lfailure ) return
c                 --- fill in the rest of the atoms using CBMC
                  iutry = 0
                  movetype = 5
               elseif ( initstyle .eq. 'unit cell' ) then
c                 --- nothing to do here, already done up above
               else
                  write(6,*) 'INITCONF: unknown initstyle:',initstyle
                  write(6,*) ' for molecule:',imolty,' in box:',ibox
                  lfailure = .true.
                  return
               endif

               if ( lgrowcbmc ) then
c                 --- get an initial structure using cbmc
c                 --- start parameters were set above
c                 --- set up the growth pattern for the rest of the
c                 --- molecule
                  iunit = nunit
                  call twh_schedule(lfailure,iunit,imolty,isteps,iutry
     &                 ,idum,0,movetype)
                  if ( lfailure ) return
c                 --- temporarily turn off the fields to generate the
c                 --- trial structure
                  call twh_nfield(GLB_GET,fieldtemp)
C                 --- store the true solvation type
                  call twh_isolvtype(GLB_GET,isolvtemp)
                  call twh_nfield(GLB_SET,0)
c                 --- reset the solvation type to none
                  call twh_isolvtype(GLB_SET,SOLV_NONE)
c                 --- perform the molecule growth
                  isafety = 0
                  lterm = .true.
                  write(6,'(a46,i5)') 
     &                 'INITCONF: building template for molecule type:'
     &                 ,imolty
c                  do while ( lterm .and. isafety .lt. 50 )
                  do while ( lterm .and. isafety .lt. 1 )
                     isafety = isafety + 1
                     write(6,'(a31,i5)') 
     &                    'Initial Structure Growth Trial:',isafety
                     call twh_resetcbmc(1,wrosen,vrostot,vrosvib,vrosang
     &                    ,vrostor,vrosext,vrosintra,vrosinter
     &                    ,vroselect,vrosewald)
                     call twh_rosenbluth(lfailure,.true.,lterm,testcord
     &                    ,CRD_REAL,1,imolty,isteps,ibox,iunit,dpscale
     &                    ,wrosen(1),vrostot(1),vrosvib(1),vrosang(1)
     &                    ,vrostor(1),vrosext(1),vrosintra(1)
     &                    ,vrosinter(1),vroselect(1),vrosewald(1))
                     if ( lfailure ) return
                  enddo
c                 --- restore the fields
                  call twh_nfield(GLB_SET,fieldtemp)
c                 --- restore the solvation type
                  call twh_isolvtype(GLB_SET,isolvtemp)

c                 --- see if we got only due to safety check
                  if ( lterm ) then
                     write(6,*) 
     &                    'unable to generate an initial conformation'
                     write(6,*) 'molecule type ',imolty
                     lfailure = .true.
                     return
                  endif

c                 --- store the template
                  do iunit = 1,nunit
                     call twh_coordinates(GLB_GET,testcord,ichain,iunit
     &                    ,xtemp,ytemp,ztemp)
                     samx(imolty,iunit) = xtemp
                     samy(imolty,iunit) = ytemp
                     samz(imolty,iunit) = ztemp
                  enddo
               endif
c              --- jump out point in case we are using previous template
 100           continue
            endif
         enddo

         ntest = 0
         iadjust = 0
         call twh_inixyz(GLB_GET,ibox,inix,iniy,iniz)
         do kc = 1, iniz
            do ic = 1, inix
               do jc = 1, iniy
                  if ( initboxtype .eq. 'unit cell' ) then
c                    --- loop through each of these unit cell images
                     istart = 1
                     iend = 0
                     do imolty = 1,nmolty
                        iend = iend + cell_initmol(imolty)
                     enddo
c                    --- initialize
                     intemp = 1
                     iadjust = 0
                     uc_count = 0
                  else
                     istart = 1
                     iend = 1
                  endif

                  do itest = istart,iend
                     ntest = ntest + 1
                     call twh_nchbox(GLB_GET,ibox,nchbox)
                     if (ntest .le. nchbox ) then
                        ichain = ichain+1
                        call twh_nboxi(GLB_SET,ichain,ibox)
                     else
                        goto 102
                     endif
                     if ( initboxtype .eq. 'unit cell' ) then
                        lfound = .false.
                        do while ( .not. lfound )
                           if ( iadjust .lt. cell_initmol(intemp) ) then
                              iadjust = iadjust + 1
                              lfound = .true.
                           else
c                             --- time to go to the next moltype
                              iadjust = 0
                              intemp = intemp + 1
                           endif
                           if ( intemp .gt. nmolty ) then
                              write(6,*) 'INITCONF: failure to find'
     &                             ,' valid moltyp for unit cell'
                              lfailure = .true.
                              return
                           endif
                        enddo
                     else
c                       --- take molecule types in order
                        do imolty = 1, nmolty
                           if ( ntest .le. mcmtma(imolty,ibox) ) then
                              intemp = imolty
                              goto 19
                           endif
                        enddo
 19                     continue
                     endif
c                    --- set chain identity information
                     call twh_moltyp(GLB_SET,ichain,intemp)
                     call twh_nunit(GLB_GET,intemp,nunit)
                     call twh_ncmt(GLB_INCR,ibox,intemp,1)
c                    --- update globalpos for the next pass
                     atomcount = atomcount + nunit
                     itemp = ichain+1
                     globalpos = 3*atomcount + 1
                     call twh_globalpos(GLB_SET,itemp,globalpos)

c                    --- get the position of all atoms for this molecule
c                    --- depending only upon the initstyle and not 
c                    --- initlattice
                     do iunit = 1,nunit
                        call twh_initstyle(GLB_GET,ibox,intemp
     &                       ,initstyle)
                        if ( initstyle .eq. 'full cbmc'
     &                       .or. initstyle .eq. 'template'
     &                       .or. initstyle .eq. 'nanotube'
     &                       .or. initstyle .eq. 'partial cbmc'
     &                       .or. initstyle .eq. 'helix cbmc'
     &                       ) then
c                          --- positions stored in our template array
                           xtemp = samx(intemp,iunit)
                           ytemp = samy(intemp,iunit)
                           ztemp = samz(intemp,iunit)
                           call twh_coordinates(GLB_SET,testcord,ichain
     &                          ,iunit,xtemp,ytemp,ztemp)
                        elseif ( initstyle .eq. 'coords' ) then
c                          --- atoms read directly from template
                           read(IFILE_COORDS,*) xtemp,ytemp,ztemp
                           call twh_coordinates(GLB_SET,testcord,ichain
     &                          ,iunit,xtemp,ytemp,ztemp)
                        elseif ( initstyle .eq. 'unit cell' ) then
                           uc_count = uc_count + 1
                           call twh_coordinates(GLB_UGT,CRD_TEMP,idum,
     &                          uc_count,xtemp,ytemp,ztemp)
c                          --- add in the extra displacement for the 
c                          --- unit cell
                           xtemp = xtemp + dble(ic-1)*cell_hmatrix(1,1)
     &                          + dble(jc-1)*cell_hmatrix(2,1)
     &                          + dble(kc-1)*cell_hmatrix(3,1)
                           ytemp = ytemp + dble(ic-1)*cell_hmatrix(1,2)
     &                          + dble(jc-1)*cell_hmatrix(2,2)
     &                          + dble(kc-1)*cell_hmatrix(3,2)
                           ztemp = ztemp + dble(ic-1)*cell_hmatrix(1,3)
     &                          + dble(jc-1)*cell_hmatrix(2,3)
     &                          + dble(kc-1)*cell_hmatrix(3,3)
c                          --- store our completely adjusted coords
                           call twh_coordinates(GLB_SET,testcord,ichain
     &                          ,iunit,xtemp,ytemp,ztemp)
                        else
                           write(6,*) 'INITCONF: unknown initstyle'
                           write(6,*) ibox,intemp,initstyle
                           lfailure = .true.
                           return
                        endif
                     enddo
c                    --- adjust the molecule coordinates according to 
c                    --- initlattice
                     call twh_initlattice(GLB_GET,ibox,intemp
     &                    ,initlattice)
                     if ( initlattice .eq. 'center' ) then
c                       --- compute the COM of molecule template
                        call twh_ctrmas(lfailure,testcord,ibox,ichain
     &                       ,CTR_INITCONF)
                        if ( lfailure ) return
c                       --- compute the center of the simulation box
                        arba = 0.5d0
                        arbb = 0.5d0
                        arbc = 0.5d0
                        call twh_arbtocart(ibox,arba,arbb,arbc,xcen,ycen
     &                       ,zcen)
c                       --- compute the offset between the center of the
c                       --- box and the center of mass of the template
                        call twh_com(GLB_GET,CRD_REAL,ichain
     &                       ,xcmi,ycmi,zcmi)
                        xoffset = xcen - xcmi
                        yoffset = ycen - ycmi
                        zoffset = zcen - zcmi
c                       --- adjust the positions
                        do iunit = 1,nunit
                           call twh_coordinates(GLB_GET,testcord,ichain
     &                          ,iunit,xtemp,ytemp,ztemp)
                           xtemp = xtemp + xoffset
                           ytemp = ytemp + yoffset
                           ztemp = ztemp + zoffset
                           call twh_coordinates(GLB_SET,testcord,ichain
     &                          ,iunit,xtemp,ytemp,ztemp)
                        enddo
                     elseif ( initlattice .eq. 'none') then
c                       --- nothing needs to be done here, do not adjust
                     elseif ( initlattice .eq. 'simple cubic' ) then
c                       --- atom coordinates evenly spaced in box
                        arba = dble(ic)*ux(ibox) - 0.5d0*ux(ibox)
                        arbb = dble(jc)*uy(ibox) - 0.5d0*uy(ibox)
                        arbc = dble(kc)*uz(ibox) - 0.5d0*uz(ibox)
c                       --- transform to cartesian coordinates
                        call twh_arbtocart(ibox,arba,arbb,arbc
     &                       ,xoffset,yoffset,zoffset)
c                       --- displace all of the atoms
                        do iunit = 1,nunit
                           call twh_coordinates(GLB_GET,testcord,ichain
     &                          ,iunit,xtemp,ytemp,ztemp)
                           xtemp = xtemp + xoffset
                           ytemp = ytemp + yoffset
                           ztemp = ztemp + zoffset
                           call twh_coordinates(GLB_SET,testcord,ichain
     &                          ,iunit,xtemp,ytemp,ztemp)
                        enddo
                     else
                        write(6,*) 'INITCONF: unknown initlattice'
                        write(6,*) 'ibox,intemp,initlattice',ibox,intemp
     &                       ,initlattice
                        lfailure = .true.
                        return
                     endif

c                    --- setup the new coordinates
                     do iunit = 1,nunit
                        call twh_coordinates(GLB_GET,testcord,ichain
     &                       ,iunit,xtemp,ytemp,ztemp)
                        call twh_coordinates(GLB_SET,CRD_REAL,ichain
     &                       ,iunit,xtemp,ytemp,ztemp)
                     enddo
                  enddo
               enddo
            enddo
         enddo
 102  continue

c     --- close the towhee_template file if it was opened
      if ( lopentemplate ) close(IFILE_TEMPLATE)

      if ( initboxtype .eq. 'unit cell' ) then
c        --- make sure that the number of molecules assigned via the 
c        --- unit cells matches the number of molecules set by nchain
         if ( ichain .ne. nchain ) then
            write(6,*) 'INITCONF: mismatch of nchain and the number'
     &           ,' of molecules implied by unit cell construction'
            write(6,*) 'nchain:',nchain
            write(6,*) 'implied nchain:',ichain
            lfailure = .true.
            return
         endif
      endif

c     --- need to setup the data structures for any molecules that 
c     --- are residing in the ideal gas resevoir
      call twh_ensemble(GLB_GET,ensemble)
      if ( ensemble .eq. ENS_muVT ) then
         ibox = 0
         idealchain = ichain
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            call twh_nmolectyp(GLB_GET,imolty,nmtype)
            do boxtest = 1,numboxes
               call twh_initmol(GLB_GET,boxtest,imolty,initmol)
               nmtype = nmtype - initmol
            enddo
            do ntest = 1,nmtype
               idealchain = idealchain + 1
               call twh_nboxi(GLB_SET,idealchain,ibox)
               call twh_moltyp(GLB_SET,idealchain,imolty)
               call twh_ncmt(GLB_INCR,ibox,imolty,1)
               call twh_nchbox(GLB_INCR,ibox,1)
c              --- update globalpos for the next pass
               atomcount = atomcount + nunit
               itemp = idealchain+1
               globalpos = 3*atomcount + 1
               call twh_globalpos(GLB_SET,itemp,globalpos)
c              --- no coordinates to set
            enddo
         enddo
      endif

c     --- close the towhee_coords file if it was opened
      if ( lopencoords ) then
         rewind(IFILE_COORDS)
         close(IFILE_COORDS)
      endif
      if ( lopenpartial ) close(IFILE_PARTIAL)

c      write(6,*) 'finished INITCONF'

      return
      end


      subroutine twh_ini_read_unitcell(lfailure,cell_hmatrix
     &     , cell_initmol)
c     ******************************************************************
c     * reads the file 'towhee_cell'.  The global variables hmatrix    *
c     * and initmol are set.  Coordinates are pushed into              *
c     * tempcoordinates stack.                                         *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  cell_hmatrix                                                  *
c     *  cell_initmol                                                  *
c     *                                                                *
c     * Error conditions:                                              *
c     *   failure to open file                                         *
c     *                                                                *
c     * split from original engmolec 02-14-2006 by MAW                 *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      double precision cell_hmatrix
      dimension cell_hmatrix(3,3)
      integer cell_initmol
      dimension cell_initmol(NTMAX)
c     --- local variables
      character*15 idname
      integer idim,jdim,idum,ibox,imolty,nmolty,iunit,numboxes
      integer uc_count,imolec,nunit,inix,iniy,iniz,initmol
      double precision xtemp,ytemp,ztemp,dvalue

c     --- retrieve constants
      call twh_numboxes(GLB_GET,numboxes)

c     --- this format is substantially different from the others
c     --- we need to open the unit cell file
      idname = 'initconf'
      call twh_nmolty(GLB_GET,nmolty)

      call twh_openfile(lfailure,IFILE_UNITCELL,'towhee_cell','DEFAULT',
     &      'OLD','DEFAULT')
      if ( lfailure ) return

      call twh_read_labeled_integer_array(lfailure,IFILE_UNITCELL,
     &   'cell_initmol',nmolty, .true., idname, cell_initmol)
      if ( lfailure ) return
      call twh_read_label(lfailure,IFILE_UNITCELL,'cell_hmatrix'
     &     ,.false.,idname)
      if ( lfailure ) return
      do idim = 1,3
         read(IFILE_UNITCELL,*) (cell_hmatrix(idim,jdim),jdim=1,3)
         write(6,'(a13,3(1x,f10.5))') 'cell_hmatrix:'
     &        ,(cell_hmatrix(idim,jdim),jdim=1,3)
      enddo

c        --- read the coordinates for the unit cell and stash them 
      call twh_read_label(lfailure,IFILE_UNITCELL,'cell_coordinates'
     &     ,.false.,idname)
      if ( lfailure ) return
      uc_count = 0
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do imolec = 1,cell_initmol(imolty)
            do iunit = 1,nunit
c                 --- increment count of atoms
               uc_count = uc_count + 1
               read(IFILE_UNITCELL,*) xtemp,ytemp,ztemp
c                 --- stash these coordinates in tempcoords
               call twh_coordinates(GLB_UST,CRD_TEMP,idum,uc_count
     &              ,xtemp,ytemp,ztemp)
            enddo
         enddo
      enddo
c        --- set all of the true arrays as multiples of the cells
      do ibox = 1,numboxes
         call twh_inixyz(GLB_GET,ibox,inix,iniy,iniz)
         do imolty = 1,nmolty
            initmol = inix*iniy*iniz*cell_initmol(imolty)
            call twh_initmol(GLB_SET,ibox,imolty,initmol)
         enddo
         do jdim = 1,3
            dvalue = inix*cell_hmatrix(1,jdim)
            call twh_hmatrix(GLB_SET,ibox,1,jdim,dvalue)
            dvalue = iniy*cell_hmatrix(2,jdim)
            call twh_hmatrix(GLB_SET,ibox,2,jdim,dvalue)
            dvalue = iniz*cell_hmatrix(3,jdim)
            call twh_hmatrix(GLB_SET,ibox,3,jdim,dvalue)
         enddo
      enddo
c        --- close this file
      close(IFILE_UNITCELL)

      return
      end
