#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2011 Marcus G. Martin                       *
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
      subroutine twh_writedynamo(lfailure,ntypes)
c     ******************************************************************
c     * output setfl and funcfl format eam functions to files          *
c     *                                                                *
c     * originally submitted 03-10-2004 by A.P. Thompson               *
c     * last modified 08-02-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EMBED
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_LHERE
#define FUNCTION_GET_MASS
#define FUNCTION_VEMBED
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      logical lovrlap,lzero
c     --- character string scalars
      character*2 ename
      character*10 eleflag
      character*20 eamembedstyle
      character*30 filename
c     --- integer scalars
      integer eam_nrho,eam_nr,ieam,iscan
      integer numone,numtwo,itype,ihere
      integer itest,itestmax,ioutput,ielement
c     --- double precision scalars
      double precision density
      double precision eam_drho,eam_dr,eam_rcutmax,eam_arbmax
      double precision rij,vpair,compare,value,rijsq,bij,bji
      double precision dummy,rcut
      double precision rtest,dv
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c     --- make sure we have an EAM forcefield
      if ( twh_get_classical_potential()
     &     .ne. 'Embedded Atom Method' ) return

c     --- retrieve constants
      call twh_rcut(GLB_GET,rcut)

c     --- hand set variables for the increment of the output
      eam_drho = 0.01
      eam_dr = 0.01

c     --- rcut is now determined in setclassical, used that here
      eam_rcutmax = rcut
      eam_nr = int(eam_rcutmax/eam_dr)

c     --- determine the largest density unit value
      eam_arbmax = 0.0d0
      do itype = 1,ntypes
         call twh_eam_embed_style(GLB_GET,itype,eamembedstyle)
         if ( eamembedstyle .eq. 'table' ) then
c           --- can find the maximum density from the final table entry
            call twh_eam_nembed(GLB_GET,itype,iscan)
            call twh_eam_embed(GLB_GET,1,itype,iscan,compare)
         else
c           --- have to do this the hard way
c           --- for now, find the maximum density with distance and 
c           --- multiply by 12
            itestmax = 100
            compare = 0.0d0
            do itest = 1,itestmax
               rtest = (dble(itest)/dble(itestmax))*rcut
c              --- ignore low distance values as they don't occur
               if ( rtest .gt. 1.0d0 ) then
                  call twh_bondorder(lfailure,itype,itype,rtest,bij,bji
     &                 ,lzero)
                  if ( lfailure ) return
                  if ( .not. lzero ) then
                     if ( bij .gt. compare ) compare = bij
                     if ( bji .gt. compare ) compare = bji
                  endif
               endif
            enddo
c           --- multiply this by 12 just to be safe
            compare = 12.0d0*compare
         endif
         if ( compare .gt. eam_arbmax ) eam_arbmax = compare
      enddo
      eam_nrho = int(eam_arbmax/eam_drho)

c     --- output the potential in the setfl format
      filename = 'setfl_out'
      call twh_openfile(lfailure,97,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return
      write(97,*) '*************************************************'
      write(97,*) 'setfl format force field file generated by TOWHEE'
      write(97,*) '*************************************************'
      write(97,*) ntypes
      write(97,*) eam_nrho,eam_drho,eam_nr,eam_dr,eam_rcutmax

      do numone = 1,ntypes
c        --- get the element number
         eleflag = 'get int'
         ename = twh_get_elementname(numone)
         call twh_getelement(lfailure,eleflag,dummy,ename,ielement)
         if ( lfailure ) return
         write(97,*) ielement, twh_get_mass(numone)

         density = 0.0
         do ieam = 1,eam_nrho
            value = twh_vembed(lfailure,numone,density)/CON_EVTOK
            if ( lfailure ) return
            write(97,'(e30.20)') value 
            density = density + eam_drho
         enddo

         rij = 0.0
         do ieam = 1,eam_nr
            call twh_bondorder(lfailure,numone,numone,rij,bij,bji,lzero)
            if ( lfailure ) return
            density = bij
            write(97,'(e30.20)') density
            rij = rij+eam_dr
         enddo
      enddo

      do numone = 1,ntypes
         do numtwo = 1,numone
            itype = twh_getnbtype(numone,numtwo)
            rij = 0.0
            do ieam = 1,eam_nr
               rijsq = rij*rij
               vpair = twh_vtwobody(lfailure,rijsq,numone,numtwo,0
     &              ,lovrlap,vf,dv)
               if ( lfailure ) return
               value = rij*vpair/CON_EVTOK
               write(97,'(e30.20)') value
               rij = rij+eam_dr
            enddo
         enddo
      enddo
      close(97)

c     --- find out how many types we have here
      ihere = 0
      write(6,*) 'ntypes',ntypes
      do numone = 1,ntypes
         write(6,*) 'numone,lhere',numone,twh_get_lhere(numone)
         if ( twh_get_lhere(numone) ) then
            ihere = ihere + 1
            ioutput = numone
         endif
      enddo

      if ( ihere .eq. 1 ) then
c        --- output the potential in the funcfl format
         filename = 'funcfl_out'
         call twh_openfile(lfailure,98,filename,'DEFAULT','UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) return
c        --- Line 1, comment 
         write(98,*) 'funcfl format file generated by Towhee'
c        --- Line 2, ielem, amass, blat, lat
         eleflag = 'get int'
         ename = twh_get_elementname(ioutput)
         call twh_getelement(lfailure,eleflag,dummy,ename,ielement)
         if ( lfailure ) return
         write(6,*) 'element,ielement',ename,ielement
         stop
         write(98,*) ''
         close(98)
      endif

      return
      end
