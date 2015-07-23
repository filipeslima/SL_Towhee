c     $Id: leastsquares.F,v 1.10 2006/02/09 20:47:36 mattwycz Exp $
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2005 Marcus G. Martin                       *
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
      subroutine twh_leastsquares(ndatapoints,nvariables,maxdp,maxv
     &     ,zzz,aaa,uuu,lscalefit,lscale)
c     ******************************************************************
c     * performs a linear leastsquares fit for the equation            *
c     * F = Sum{ndatapoints} ( Sum{nvariables} aaa*zzz - uuu )^2       *
c     *                                                                *
c     * originally written 08-15-2003 by M.G. Martin                   *
c     * last modified 10-18-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ndatapoints,nvariables,maxdp,maxv
      double precision aaa,zzz,uuu
      logical lscalefit,lscale
      dimension aaa(maxv)
      dimension uuu(maxdp)
      dimension zzz(maxv,maxdp)
      dimension lscale(maxdp)
c     ---- local variables
      double precision zwork,invzwork
      dimension zwork(nvariables,nvariables)
      dimension invzwork(nvariables,nvariables)
      double precision uwork
      dimension uwork(nvariables)
      integer idex,jdex,idp,iidex,jjdex,nscalevariables
      logical lmaterror

c     --- initialize arrays
      do idex = 1,nvariables
         aaa(idex) = 0.0d0
         uwork(idex) = 0.0d0
         do jdex = 1,nvariables
            zwork(idex,jdex) = 0.0d0
         enddo
      enddo

c     --- set up the uwork array = Sum z(i)*u
      if (lscalefit) then
         nscalevariables = 0
         jdex = 0
         do idex = 1,nvariables
            if (lscale(idex)) then
               jdex = jdex + 1
               nscalevariables = nscalevariables + 1
               do idp = 1,ndatapoints
                  uwork(jdex) = uwork(jdex) + zzz(jdex,idp)*uuu(idp)
               enddo
            endif
         enddo
      else
         do idex = 1,nvariables
            do idp = 1,ndatapoints
               uwork(idex) = uwork(idex) + zzz(idex,idp)*uuu(idp)
            enddo
         enddo
      endif

c     --- set up the zwork array = Sum z(i)*z(j)
      if (lscalefit) then
         iidex = 0
         do idex = 1,nvariables   
            if (lscale(idex)) then
               iidex = iidex + 1
               jjdex = 0
               do jdex = 1,nvariables
                  if (lscale(jdex)) then
                     jjdex = jjdex + 1
                     do idp = 1,ndatapoints
                        zwork(iidex,jjdex) = zwork(iidex,jjdex) 
     &                     + zzz(iidex,idp)*zzz(jjdex,idp)
                     enddo
                  endif
               enddo
            endif
         enddo
      else
         do idex = 1,nvariables   
            do jdex = 1,nvariables
               do idp = 1,ndatapoints
                  zwork(idex,jdex) = zwork(idex,jdex) 
     &                     + zzz(idex,idp)*zzz(jdex,idp)
               enddo
            enddo
         enddo
      endif

      if (.not. lscalefit) nscalevariables = nvariables
c     --- invert zwork
      call twh_invert(nvariables,nscalevariables,zwork,invzwork
     &     ,lmaterror)

      if ( lmaterror ) then
         do idex = 1,nvariables
            write(6,*) 'idex,zwork',idex
     &           ,(zwork(idex,jdex),jdex=1,nvariables)
         enddo
      endif

c     --- compute aaa
      do idex = 1,nscalevariables
         do jdex = 1,nscalevariables
            aaa(idex) = aaa(idex) + uwork(jdex)*invzwork(jdex,idex)
         enddo
      enddo

      return
      end


      subroutine twh_invert( n, nscale, A, Ainv, lmaterror )
c     ******************************************************************
c     * routine to invert a matrix by Gaussian elimination             *
c     * Ainv=inverse(A)    A(n,n)   Ainv(n,n)                          *
c     *                                                                *
c     * originally rewritten from web version 08-15-2003 by M.G. Martin*
c     * (http://www.nauticom.net/www/jdtaft/FortranMatrix.htm)         *
c     * last modified 10-18-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lmaterror
      integer n,nscale
      double precision A,Ainv
      dimension A(n,n),Ainv(n,n)
c     --- local variables
      integer n2,i,j,k
      double precision aaa,bbb
      double precision ddd
      dimension ddd(nscale,2*nscale)
c     --- initialize error condition
      lmaterror = .false.

c     --- initialize the reduction matrix
      n2 = 2*nscale
      do i = 1,nscale
         do j = 1,nscale
            ddd(i,j) = A(i,j)
            ddd(i,nscale+j) = 0.
         enddo
         ddd(i,nscale+i) = 1.
      enddo
      
c     --- do the reduction 
      do i = 1,nscale
         aaa = ddd(i,i)
         if (aaa .eq. 0.) then
            write(6,*) '*** ERROR: Singular matrix ***'
            write(6,*) 'Zero value diagonal',i,aaa
            lmaterror = .true.
            return
         endif
         do j = 1,n2
            ddd(i,j) = ddd(i,j)/aaa
         enddo
         do k = 1,nscale
            if ( (k-i) .ne. 0 ) then
               bbb = ddd(k,i)
               do j = 1,n2
                  ddd(k,j) = ddd(k,j) - bbb*ddd(i,j)
               enddo
            endif
         enddo
      enddo
      
c     --- copy result into output matrix
      do i = 1,nscale
         do j = 1,nscale
            Ainv(i,j) = ddd(i,j+nscale)
         enddo
      enddo
      
      return
      end

