#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2013 Marcus G. Martin                       *
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
      function twh_vembed(lfailure,ntii,density)
c     ******************************************************************
c     * computes the embedding function for the EAM potential          *
c     *                                                                *
c     * originally written 11-19-2003 by  M.G. Martin                  *
c     * last update 07-22-2013 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_LT
#define FUNCTION_VEMBED
#define FUNCTION_EXTRACTEMBED
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure
      integer ntii
      double precision density
c     --- local variables
      integer eamindex,eamembednum
      double precision dsq,embedone,embedtwo,embedthr
      double precision rhomin,rhozero,eamalpha,eambeta,densdiff
      double precision belm,beln,bela,belb,belc,densratio
c     --- double precision arrays
      double precision belrho
      dimension belrho(0:9)

      call twh_eam_embed_num(GLB_GET,ntii,eamembednum)
      if ( eamembednum .eq. EAM_EMBED_BELASH ) then
         call twh_eam_embed(GLB_GET,1,ntii,5,rhomin)
         if ( density .lt. rhomin ) then
c           --- if rho < rhomin then v = alpha*sqrt(rho) + beta*rho
            call twh_eam_embed(GLB_GET,1,ntii,4,eamalpha)
            call twh_eam_embed(GLB_GET,2,ntii,4,eambeta)
            twh_vembed = eamalpha*dsqrt(density) + eambeta*density
         else
c           --- a0 + a1*(rho-rho0)^2 + a3*(rho-rho0)^3
            call twh_eam_embed(GLB_GET,2,ntii,5,rhozero)
            call twh_eam_embed(GLB_GET,2,ntii,1,embedone)
            call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
            call twh_eam_embed(GLB_GET,2,ntii,3,embedthr)
            densdiff = density - rhozero
            dsq = densdiff*densdiff
            twh_vembed = embedone + embedtwo*(dsq)
     &           + embedthr*dsq*densdiff
         endif
      elseif ( eamembednum .eq. EAM_EMBED_BELASHTRHO ) then
c        --- Belaschenko piecewise potential based on the density
         call twh_eam_embed(GLB_GET,1,ntii,11,belrho(0))
         call twh_eam_embed(GLB_GET,2,ntii,11,belrho(1))
         call twh_eam_embed(GLB_GET,3,ntii,11,belrho(2))
         call twh_eam_embed(GLB_GET,1,ntii,12,belrho(3))
         call twh_eam_embed(GLB_GET,2,ntii,12,belrho(4))
         call twh_eam_embed(GLB_GET,3,ntii,12,belrho(5))
         call twh_eam_embed(GLB_GET,1,ntii,13,belrho(6))
         call twh_eam_embed(GLB_GET,2,ntii,13,belrho(7))
         call twh_eam_embed(GLB_GET,3,ntii,13,belrho(8))
         call twh_eam_embed(GLB_GET,1,ntii,14,belrho(9))
c        --- go small to large density
         if ( twh_cmp_lt(density,belrho(7)) ) then
c           --- 0 to belrho(7)
            call twh_eam_embed(GLB_GET,1,ntii,8,bela)
            call twh_eam_embed(GLB_GET,2,ntii,8,belb)
            call twh_eam_embed(GLB_GET,3,ntii,8,belc)
            densdiff = density-belrho(7)
            densratio = density/belrho(7)
            twh_vembed = (bela + belb*densdiff + belc*densdiff*densdiff)
     &           *(2.0d0*densratio - densratio*densratio)
         elseif ( twh_cmp_lt(density,belrho(1)) ) then
c           --- somewhere in belrho(7) to belrho(1)
c           --- this group uses the same function, find the subinterval
            if ( twh_cmp_lt(density,belrho(6)) ) then
               eamindex = 7
            elseif ( twh_cmp_lt(density,belrho(5)) ) then
               eamindex = 6
            elseif ( twh_cmp_lt(density,belrho(4)) ) then
               eamindex = 5
            elseif ( twh_cmp_lt(density,belrho(3)) ) then
               eamindex = 4
            elseif ( twh_cmp_lt(density,belrho(2)) ) then
               eamindex = 3
            else
               eamindex = 2
            endif
            densdiff = density - belrho(eamindex-1)
            call twh_eam_embed(GLB_GET,1,ntii,eamindex,bela)
            call twh_eam_embed(GLB_GET,2,ntii,eamindex,belb)
            call twh_eam_embed(GLB_GET,3,ntii,eamindex,belc)
            twh_vembed = bela + densdiff*(belb + belc*densdiff)
         elseif ( twh_cmp_lt(density,belrho(8)) ) then
c           --- this interval doesn't need a b and refers to another
c           --- central rho that is between belrho(1) and belrho(8)
            call twh_eam_embed(GLB_GET,1,ntii,1,bela)
            call twh_eam_embed(GLB_GET,3,ntii,1,belc)
            densdiff = density - belrho(0)
            twh_vembed = bela + belc*densdiff*densdiff
         elseif ( twh_cmp_lt(density,belrho(9)) ) then
c           --- between belrho(8) and belrho(9), getting to high density
            call twh_eam_embed(GLB_GET,1,ntii,9,bela)
            call twh_eam_embed(GLB_GET,2,ntii,9,belb)
            call twh_eam_embed(GLB_GET,3,ntii,9,belc)
            call twh_eam_embed(GLB_GET,2,ntii,14,belm)
            densdiff = density - belrho(8)
            twh_vembed = bela + belb*densdiff + belc*(densdiff**belm)
         else
c           --- very high density, beyond belrho(9)
            call twh_eam_embed(GLB_GET,1,ntii,10,bela)
            call twh_eam_embed(GLB_GET,2,ntii,10,belb)
            call twh_eam_embed(GLB_GET,3,ntii,10,belc)
            call twh_eam_embed(GLB_GET,3,ntii,14,beln)
            densdiff = density - belrho(9)
            twh_vembed = bela + belb*densdiff + belc*(densdiff**beln)
         endif
      elseif ( eamembednum .eq. EAM_EMBED_LOGARITHMIC ) then
c        --- functional form from JJ Hoyt and J Zimmerman
c        --- f(density) = density*(a*ln(density)+b)
         call twh_eam_embed(GLB_GET,2,ntii,1,embedone)
         call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
         twh_vembed = density*((embedone*log(density))+embedtwo)
      elseif ( eamembednum .eq. EAM_EMBED_POWHT ) then
c        --- Ackland two terms of square root and square
         call twh_eam_embed(GLB_GET,2,ntii,1,embedone)
         call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
         twh_vembed = embedone * dsqrt(density) 
     &        + embedtwo*density*density 
      elseif ( eamembednum .eq. EAM_EMBED_POWHTF ) then
c        --- power series with 0.5, 2, and 4 powers
         dsq = density*density
         call twh_eam_embed(GLB_GET,2,ntii,1,embedone)
         call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
         call twh_eam_embed(GLB_GET,2,ntii,3,embedthr)
         twh_vembed = embedone * dsqrt(density)
     &        + embedtwo*dsq + embedthr*dsq*dsq
      elseif ( eamembednum .eq. EAM_EMBED_SR ) then
c        --- functional form from JJ Hoyt and J Zimmerman
c        --- f(density) = b*density - a*(sqrt(density))
         call twh_eam_embed(GLB_GET,2,ntii,1,embedone)
         call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
         twh_vembed = (embedtwo*density) - (embedone*dsqrt(density))
      elseif ( eamembednum .eq. EAM_EMBED_TABLE ) then
         twh_vembed = twh_extractembed(lfailure,INTERP_NORMAL
     &        ,ntii,density)
         if ( lfailure ) return
      else
         write(6,*) 'VEMBED: invalid eam_embed_style'
         write(6,*) 'ntii',ntii
     &        ,'eam_embed_num:',eamembednum
         stop
      endif

c      write(6,*) 'finished VEMBED'

      return
      end
