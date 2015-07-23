c     --- toggles for the global variables
#define GLB_SET 2001
#define GLB_GET 2002
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2008 Marcus G. Martin                       *
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
      function dranlux()
c     ******************************************************************
c     * this function has been rewritten to generate pseudorandom      *
c     * numbers just for the RANLUX algorithm using a modified version *
c     * of the ranlux generator as posted on the web by Phillip Helbiq *
c     * Note: we use the luxury levels defined in                      *
c     * F. James; "RANLUX: A Fortran implementation of the high-quality*
c     * pseudorandom number generator of Luscher", Comp. Phys. Comm.   *
c     * 79 111-114 (1994).                                             *
c     * where he is utilizing the algorithm of                         *
c     * M. Luscher; "A portable high-quality random number generator   *
c     * for lattice field-theory simulations" Comp. Phys. Comm.        *
c     * 79 100-110 (1994).                                             *
c     *                                                                *
c     * completely rewritten 08-07-2003 by M.G. Martin                 *
c     * last modified 08-16-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- function declaration
c     --- this is an exception from the normal Towhee format of 
c     --- declaring all functions in functions.h do avoid having any of 
c     --- the include files in this routine in case they would conflict
c     --- with the ranlux implementation
c     --- Note: random is including in the function list in functions.h 
c     --- for compatability with all of the other routines
      double precision dranlux
      real*8 stuff
      dimension stuff(1)

      call ranlux(stuff,1)
      dranlux = stuff(1)

      return
      end

C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C          
C       references: 
C  M. Luscher, Computer Physics Communications  79 (1994) 100
C  F. James, Computer Physics Communications 79 (1994) 111
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,k2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and k2   ++
C!!!               should be set to zero unless restarting at a break++ 
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,k2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and k2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+k2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     * last modified 08-02-2011 by M.G. Martin                       *
c     *****************************************************************
      subroutine ranlux(rvec,lenv)
      implicit none
c     --- definitions
#define LXDFLT 3
#define MAXLEV 4
#define TWOP12 4096.
#define IGIGA 1000000000
#define JSDFLT 314159265
#define ITWO24 2**24
#define ICONS 2147483563
c     --- variables passed to/from the subroutine
      integer lenv
      real*8 rvec
      dimension rvec(lenv)
c     --- local variables
C     default
C     Luxury Level   0     1     2   *3*    4
C     corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C     1    1.5    2     3     5   on fast mainframe
C     
c     --- logical scalars
      logical notyet
c     --- 32-bit integer scalars
      integer*4 inseed,nskip,luxlev,ranlux_ndskip,intwentyfour,mkount
      integer*4 kount,jseed,lp,i,k,itwentyfour,jtwentyfour,ivec,isk
c     --- 32-bit integer arrays
      integer*4 next
      dimension next(24)
      integer*4 iseeds
      dimension iseeds(24)
c     --- real scalars
      real*8 carry,twomtwentyfour,twomtwelve,uni
c     --- real arrays
      real*8 seeds
      dimension seeds(24)

      call ranlux_notyet(GLB_GET,notyet)
C     --- NOTYET is .TRUE. if no initialization has been performed yet.
C     --- Default Initialization by Multiplicative Congruential
      if ( notyet ) then
         notyet = .FALSE.
         call ranlux_notyet(GLB_SET,notyet)
         jseed = JSDFLT  
         inseed = jseed
         call ranlux_inseed(GLB_SET,inseed)
         write(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',jseed
         luxlev = LXDFLT
         call ranlux_luxlev(GLB_SET,luxlev)
         nskip = ranlux_ndskip(luxlev)
         call ranlux_nskip(GLB_SET,nskip)
         lp = nskip + 24
         intwentyfour = 0
         call ranlux_intwentyfour(GLB_SET,intwentyfour)
         kount = 0
         call ranlux_kount(GLB_SET,kount)
         mkount = 0
         call ranlux_mkount(GLB_SET,mkount)
         write(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
     +        luxlev,'      p =',lp
         twomtwentyfour = 1.
         do I= 1, 24
            twomtwentyfour = twomtwentyfour * 0.5
            K = jseed/53668
            jseed = 40014*(jseed-K*53668) -K*12211
            if (jseed .LT. 0)  jseed = jseed+ICONS
            iseeds(I) = MOD(jseed,ITWO24)
         enddo
         twomtwelve = twomtwentyfour * 4096.
         do I= 1,24
            seeds(I) = REAL(iseeds(I))*twomtwentyfour
            call ranlux_seeds(GLB_SET,i,seeds(i))
            next(I) = I-1
            call ranlux_next(GLB_SET,i,next(i))
         enddo
         next(1) = 24
         call ranlux_next(GLB_SET,1,next(1))
         itwentyfour = 24
         call ranlux_itwentyfour(GLB_SET,itwentyfour)
         jtwentyfour = 10
         call ranlux_jtwentyfour(GLB_SET,jtwentyfour)
         carry = 0.
         if (SEEDS(24) .EQ. 0.) carry = twomtwentyfour
         call ranlux_carry(GLB_SET,carry)
         call ranlux_twomtwentyfour(GLB_SET,twomtwentyfour)
         call ranlux_twomtwelve(GLB_SET,twomtwelve)
      endif
C     
C     The Generator proper: "Subtract-with-borrow",
C     as proposed by Marsaglia and Zaman,
C     Florida State University, March, 1989
C     

      call ranlux_nskip(GLB_GET,nskip)
      call ranlux_itwentyfour(GLB_GET,itwentyfour)
      call ranlux_jtwentyfour(GLB_GET,jtwentyfour)
      call ranlux_carry(GLB_GET,carry)
      call ranlux_twomtwentyfour(GLB_GET,twomtwentyfour)
      call ranlux_twomtwelve(GLB_GET,twomtwelve)
      call ranlux_intwentyfour(GLB_GET,intwentyfour)
      call ranlux_kount(GLB_GET,kount)

      do ivec= 1, lenv
         call ranlux_seeds(GLB_GET,itwentyfour,seeds(itwentyfour))
         call ranlux_seeds(GLB_GET,jtwentyfour,seeds(jtwentyfour))
         uni = seeds(jtwentyfour) - seeds(itwentyfour) - carry 
         if (uni .LT. 0.)  then
            uni = uni + 1.0
            carry = twomtwentyfour
         else
            carry = 0.
         endif
         seeds(itwentyfour) = uni
         call ranlux_seeds(GLB_SET,itwentyfour,seeds(itwentyfour))
         call ranlux_next(GLB_GET,itwentyfour,next(itwentyfour))
         itwentyfour = next(itwentyfour)
         call ranlux_next(GLB_GET,jtwentyfour,next(jtwentyfour))
         jtwentyfour = next(jtwentyfour)
         rvec(ivec) = uni
C        --- small numbers (with less than 12 "significant" bits) are
c        --- "padded".
         if (uni .LT. twomtwelve)  then
            call ranlux_seeds(GLB_GET,jtwentyfour,seeds(jtwentyfour))
            rvec(ivec) = rvec(ivec) + twomtwentyfour*seeds(jtwentyfour)
C           --- zero is forbidden in case someone takes a logarithm
            if (rvec(ivec) .EQ. 0.)  rvec(ivec)
     &           = twomtwentyfour*twomtwentyfour
         endif
C        --- Skipping to luxury.  As proposed by Martin Luscher.
         intwentyfour = intwentyfour + 1
         if (intwentyfour .EQ. 24)  then
            intwentyfour = 0
            kount = kount + nskip
            do isk= 1, nskip
               call ranlux_seeds(GLB_GET,itwentyfour,seeds(itwentyfour))
               call ranlux_seeds(GLB_GET,jtwentyfour,seeds(jtwentyfour))
               uni = seeds(jtwentyfour) - seeds(itwentyfour) - carry
               if (uni .LT. 0.)  then
                  uni = uni + 1.0
                  carry = twomtwentyfour
               else
                  carry = 0.
               endif
               seeds(itwentyfour) = uni
               call ranlux_seeds(GLB_SET,itwentyfour,seeds(itwentyfour))
               call ranlux_next(GLB_GET,itwentyfour,next(itwentyfour))
               itwentyfour = next(itwentyfour)
               call ranlux_next(GLB_GET,jtwentyfour,next(jtwentyfour))
               jtwentyfour = next(jtwentyfour)
            enddo
         endif
      enddo
      kount = kount + lenv
      call ranlux_mkount(GLB_GET,mkount)
      if (kount .GE. IGIGA)  then
         mkount = mkount + 1
         kount = kount - IGIGA
      endif
c     --- save new values of the parameters
      call ranlux_mkount(GLB_SET,mkount)
      call ranlux_kount(GLB_SET,kount)
      call ranlux_itwentyfour(GLB_SET,itwentyfour)
      call ranlux_jtwentyfour(GLB_SET,jtwentyfour)
      call ranlux_carry(GLB_SET,carry)
      call ranlux_twomtwentyfour(GLB_SET,twomtwentyfour)
      call ranlux_intwentyfour(GLB_SET,intwentyfour)
      return
      end

      subroutine RLUXIN(isdext)
c     ******************************************************************
c     * routine to input and float integer seeds from previous run     *
c     * if block added by Phillip Helbig after correpondence with James*
c     *                                                                *
c     * last updated 08-15-2011 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer*4 isdext
      dimension isdext(25)
c     --- local variables
c     --- logical scalars
      logical notyet
c     --- 32-bit integer scalars
      integer*4 i,isd
      integer*4 itwentyfour,jtwentyfour,intwentyfour,luxlev
      integer*4 nskip,ranlux_ndskip,inseed
c     --- 32-bit integer arrays
      integer*4 next
      dimension next(24)
c     --- real scalars
      real*8 carry,twomtwentyfour,twomtwelve
c     --- real arrays
      real*8 seeds
      dimension seeds(24)

      call ranlux_notyet(GLB_GET,notyet)
      if (notyet) then
         notyet = .FALSE.
         call ranlux_notyet(GLB_SET,notyet)
      endif
      twomtwentyfour = 1.
      do I= 1, 24
         next(I) = I-1
         call ranlux_next(GLB_SET,i,next(i))
         twomtwentyfour = twomtwentyfour * 0.5
      enddo
      next(1) = 24
      call ranlux_next(GLB_SET,1,next(1))
      twomtwelve = twomtwentyfour * 4096.
      write(6,'(A)') 
     &     'Full initialization of RANLUX with 25 Integers:'
      write(6,'(5X,5I12)') isdext
      do I= 1, 24
         seeds(I) = REAL(isdext(I))*twomtwentyfour
         call ranlux_seeds(GLB_SET,i,seeds(i))
      enddo
      carry = 0.
      if (isdext(25) .LT. 0)  carry = twomtwentyfour
      call ranlux_carry(GLB_SET,carry)
      isd = IABS(isdext(25))
      itwentyfour = MOD(isd,100)
      call ranlux_itwentyfour(GLB_SET,itwentyfour)
      isd = isd/100
      jtwentyfour = MOD(isd,100)
      call ranlux_jtwentyfour(GLB_SET,jtwentyfour)
      isd = isd/100
      intwentyfour = MOD(isd,100)
      call ranlux_intwentyfour(GLB_SET,intwentyfour)
      isd = isd/100
      luxlev = isd
      call ranlux_luxlev(GLB_SET,luxlev)
      if (luxlev .LE. MAXLEV) then
         nskip = ranlux_ndskip(luxlev)
         call ranlux_nskip(GLB_SET,nskip)
         write (6,'(A,I2)') 
     &        ' RANLUX luxury level set by RLUXIN to: ',
     +        luxlev
      elseif (luxlev .GE. 24) then
         nskip = luxlev - 24
         call ranlux_nskip(GLB_SET,nskip)
         write (6,'(A,I5)') 
     &        ' RANLUX P-VALUE SET BY RLUXIN TO:',luxlev
      else
         nskip = ranlux_ndskip(MAXLEV)
         call ranlux_nskip(GLB_SET,nskip)
         write (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',luxlev
         luxlev = MAXLEV
         call ranlux_luxlev(GLB_SET,luxlev)
      endif
      inseed = -1
      call ranlux_inseed(GLB_SET,inseed)
      call ranlux_twomtwentyfour(GLB_SET,twomtwentyfour)
      call ranlux_twomtwelve(GLB_SET,twomtwelve)

      return
      end

      subroutine RLUXUT(isdext)
c     ******************************************************************
c     * routine to ouput seeds as integers                             *
c     *                                                                *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer*4 isdext
      dimension isdext(25)
c     --- local variables
c     32-bit integer scalars
      integer*4 i,luxlev,itwentyfour,jtwentyfour,intwentyfour
c     --- real scalars
      real*8 carry
c     --- real arrays
      real*8 seeds
      dimension seeds(24)

      do I= 1, 24
         call ranlux_seeds(GLB_GET,i,seeds(i))
         isdext(I) = INT(seeds(I)*TWOP12*TWOP12)
      enddo
      call ranlux_luxlev(GLB_GET,luxlev)
      call ranlux_itwentyfour(GLB_GET,itwentyfour)
      call ranlux_jtwentyfour(GLB_GET,jtwentyfour)
      call ranlux_intwentyfour(GLB_GET,intwentyfour)
      isdext(25) = itwentyfour + 100*jtwentyfour + 10000*intwentyfour
     &     + 1000000*luxlev
      call ranlux_carry(GLB_GET,carry)
      if (carry .GT. 0.)  isdext(25) = -isdext(25)
      return
      end

      subroutine RLUXAT(lout,inout,k1,k2)
c     ******************************************************************
c     * routine to output the "convenient" restart point               *
c     *                                                                *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer*4 lout,inout,k1,k2
c     --- local variables
      integer*4 luxlev,inseed,kount,mkount

      call ranlux_luxlev(GLB_GET,luxlev)
      lout = luxlev
      call ranlux_inseed(GLB_GET,inseed)
      inout = inseed
      call ranlux_kount(GLB_GET,kount)
      k1 = kount
      call ranlux_mkount(GLB_GET,mkount)
      k2 = mkount
      return
      end

      subroutine RLUXGO(lux,ins,k1,k2)
c     ******************************************************************
c     * routine to initialize from one or three integers               *
c     *                                                                *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer*4 lux,ins,k1,k2
c     --- local variables
c     --- logical scalars
      logical notyet
c     --- 32-bit integer scalars
      integer*4 nskip,ranlux_ndskip,luxlev,intwentyfour,jseed,inseed
      integer*4 iouter,inner,izip,iziptwo,ilx,i,k,kount,mkount,isk
      integer*4 itwentyfour,jtwentyfour
c     --- 32-bit integer arrays
      integer*4 next
      dimension next(24)
      integer*4 iseeds
      dimension iseeds(24)
c     --- real scalars
      real*8 twomtwentyfour,twomtwelve,carry,uni
c     --- real arrays
      real*8 seeds
      dimension seeds(24)

      if (lux .LT. 0) then
         luxlev = LXDFLT
         call ranlux_luxlev(GLB_SET,luxlev)
      elseif (lux .LE. MAXLEV) then
         luxlev = lux
         call ranlux_luxlev(GLB_SET,luxlev)
      elseif (lux .LT. 24 .OR. lux .GT. 2000) then
         luxlev = MAXLEV
         call ranlux_luxlev(GLB_SET,luxlev)
         write (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',lux
      else
         luxlev = lux
         do ilx= 0, MAXLEV
            nskip = ranlux_ndskip(ilx)
            if (lux .EQ. nskip+24)  luxlev = ilx
         enddo
         call ranlux_luxlev(GLB_SET,luxlev)
      endif
      if (luxlev .LE. MAXLEV)  then
         nskip = ranlux_ndskip(luxlev)
         call ranlux_nskip(GLB_SET,nskip)
         write(6,'(A,I2,A,I4)') 
     &        ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
     +        luxlev,'     P=', nskip+24
      else
         nskip = luxlev - 24
         call ranlux_nskip(GLB_SET,nskip)
         write (6,'(A,I5)') 
     &        ' RANLUX P-VALUE SET BY RLUXGO TO:',luxlev
      endif
      intwentyfour = 0
      if (ins .LT. 0)  write (6,'(A)')   
     +     ' Illegal initialization by RLUXGO, negative input seed'
      if (ins .GT. 0)  then
         jseed = ins
         write(6,'(A,3I12)') 
     &        ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
     +        jseed, k1,k2
      else
         jseed = JSDFLT
         write(6,'(A)')
     &        ' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      endif
      inseed = jseed
      call ranlux_inseed(GLB_SET,inseed)
      notyet = .FALSE.
      call ranlux_notyet(GLB_SET,notyet)
      twomtwentyfour = 1.
      do I= 1, 24
         twomtwentyfour = twomtwentyfour * 0.5
         K = jseed/53668
         jseed = 40014*(jseed-K*53668) -K*12211
         if (jseed .LT. 0)  jseed = jseed+ICONS
         iseeds(I) = MOD(jseed,ITWO24)
      enddo
      twomtwelve = twomtwentyfour * 4096.
      do I= 1,24
         seeds(I) = REAL(iseeds(I))*twomtwentyfour
         call ranlux_seeds(GLB_SET,i,seeds(i))
         next(I) = I-1
         call ranlux_next(GLB_SET,i,next(i))
      enddo
      next(1) = 24
      call ranlux_next(GLB_SET,1,next(1))
      itwentyfour = 24
      call ranlux_itwentyfour(GLB_SET,itwentyfour)
      jtwentyfour = 10
      call ranlux_jtwentyfour(GLB_SET,jtwentyfour)
      carry = 0.
      if (seeds(24) .EQ. 0.) carry = twomtwentyfour
      call ranlux_carry(GLB_SET,carry)
C     If restarting at a break point, skip k1 + IGIGA*k2
C     Note that this is the number of numbers delivered to
C     the user PLUS the number skipped (if luxury .GT. 0).
      kount = k1
      mkount = k2
      if (k1+k2 .NE. 0)  then
         do iouter= 1, k2+1
            inner = IGIGA
            if (iouter .EQ. k2+1)  inner = k1
            do isk= 1, inner
               uni = seeds(jtwentyfour) - seeds(itwentyfour) - carry 
               if (uni .LT. 0.)  then
                  uni = uni + 1.0
                  carry = twomtwentyfour
               else
                  carry = 0.
               endif
               seeds(itwentyfour) = uni
               call ranlux_seeds(GLB_SET,itwentyfour,seeds(itwentyfour))
               call ranlux_next(GLB_GET,itwentyfour,next(itwentyfour))
               itwentyfour = next(itwentyfour)
               call ranlux_next(GLB_GET,jtwentyfour,next(jtwentyfour))
               jtwentyfour = next(jtwentyfour)
            enddo
         enddo
C        --- Get the right value of intwentyfour by direct calculation
         intwentyfour = MOD(kount, nskip+24)
         if (mkount .GT. 0)  then
            izip = MOD(IGIGA, nskip+24)
            iziptwo = mkount*izip + intwentyfour
            intwentyfour = MOD(iziptwo, nskip+24)
         endif
C        --- Now intwentyfour had better be between zero and 23
c        --- inclusive
         if (intwentyfour .GT. 23) then
            write (6,'(A/A,3I11,A,I5)')  
     +           '  Error in RESTARTING with RLUXGO:'
     &           ,'  The values', ins,
     +           k1, k2, ' cannot occur at luxury level', luxlev
            intwentyfour = 0
         endif
      endif
c     --- store final values
      call ranlux_itwentyfour(GLB_SET,itwentyfour)
      call ranlux_jtwentyfour(GLB_SET,jtwentyfour)
      call ranlux_carry(GLB_SET,carry)
      call ranlux_twomtwentyfour(GLB_SET,twomtwentyfour)
      call ranlux_twomtwelve(GLB_SET,twomtwelve)
      call ranlux_intwentyfour(GLB_SET,intwentyfour)
      call ranlux_mkount(GLB_SET,mkount)
      call ranlux_kount(GLB_SET,kount)
      return
      end
      
      subroutine ranlux_carry(flag,rvalue)
c     ******************************************************************
c     * the carry global variable in the random number generator       *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag
      real*8 rvalue
c     --- information variables
      real*8 carry
      save carry

      if ( flag .eq. GLB_SET ) then
         carry = rvalue
      elseif ( flag .eq. GLB_GET ) then
         rvalue = carry
      else
         write(6,*) 'RANLUX_carry: invalid global flag',flag
         stop
      endif

      return
      end
      
      subroutine ranlux_inseed(flag,ivalue)
c     ******************************************************************
c     * the inseed global variable in the random number generator      *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,ivalue
c     --- information variables
      integer*4 inseed
      save inseed

      if ( flag .eq. GLB_SET ) then
         inseed = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = inseed
      else
         write(6,*) 'RANLUX_inseed: invalid global flag',flag
         stop
      endif

      return
      end
      
      subroutine ranlux_intwentyfour(flag,ivalue)
c     ******************************************************************
c     * the intwentyfour global variable in the random number generator*
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,ivalue
c     --- information variables
      integer*4 intwentyfour
      save intwentyfour

      if ( flag .eq. GLB_SET ) then
         intwentyfour = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = intwentyfour
      else
         write(6,*) 'RANLUX_intwentyfour: invalid global flag',flag
         stop
      endif

      return
      end

      
      subroutine ranlux_itwentyfour(flag,ivalue)
c     ******************************************************************
c     * the itwentyfour global variable in the random number generator *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,ivalue
c     --- information variables
      integer*4 itwentyfour
      save itwentyfour

      if ( flag .eq. GLB_SET ) then
         itwentyfour = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = itwentyfour
      else
         write(6,*) 'RANLUX_itwentyfour: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine ranlux_jtwentyfour(flag,ivalue)
c     ******************************************************************
c     * the jtwentyfour global variable in the random number generator *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,ivalue
c     --- information variables
      integer*4 jtwentyfour
      save jtwentyfour

      if ( flag .eq. GLB_SET ) then
         jtwentyfour = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = jtwentyfour
      else
         write(6,*) 'RANLUX_jtwentyfour: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine ranlux_kount(flag,ivalue)
c     ******************************************************************
c     * the kount global variable in the random number generator       *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,ivalue
c     --- information variables
      integer*4 kount
      save kount

      if ( flag .eq. GLB_SET ) then
         kount = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = kount
      else
         write(6,*) 'RANLUX_kount: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine ranlux_luxlev(flag,ivalue)
c     ******************************************************************
c     * the luxlev global variable in the random number generator      *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,ivalue
c     --- information variables
      integer*4 luxlev
      save luxlev

      if ( flag .eq. GLB_SET ) then
         luxlev = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = luxlev
      else
         write(6,*) 'RANLUX_luxlev: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine ranlux_mkount(flag,ivalue)
c     ******************************************************************
c     * the mkount global variable in the random number generator      *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,ivalue
c     --- information variables
      integer*4 mkount
      save mkount

      if ( flag .eq. GLB_SET ) then
         mkount = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = mkount
      else
         write(6,*) 'RANLUX_mkount: invalid global flag',flag
         stop
      endif

      return
      end

      integer*4 function ranlux_ndskip(index)
c     ******************************************************************
c     * the ndskip settings for the random number generator            *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this function
      integer*4 index
c     --- 32-bit integer arrays
      integer*4 ndskip
      dimension ndskip(0:MAXLEV)
      data ndskip /0,   24,   73,  199,  365 /

      ranlux_ndskip = ndskip(index)
      return
      end

      subroutine ranlux_next(flag,index,ivalue)
c     ******************************************************************
c     * the next global variable in the random number generator        *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,index,ivalue
c     --- information variables
      integer*4 next
      dimension next(24)
      save next

      if ( flag .eq. GLB_SET ) then
         next(index) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = next(index)
      else
         write(6,*) 'RANLUX_next: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine ranlux_notyet(flag,lvalue)
c     ******************************************************************
c     * the notyet global variable in the random number generator      *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer*4 flag
      logical lvalue
c     --- information variables
      logical notyet
      save notyet
      data notyet /.true./

      if ( flag .eq. GLB_SET ) then
         notyet = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = notyet
      else
         write(6,*) 'RANLUX_notyet: invalid global flag',flag
         stop
      endif
      return
      end

      subroutine ranlux_nskip(flag,ivalue)
c     ******************************************************************
c     * the nskip global variable in the random number generator       *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,ivalue
c     --- information variables
      integer*4 nskip
      save nskip

      if ( flag .eq. GLB_SET ) then
         nskip = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = nskip
      else
         write(6,*) 'RANLUX_nskip: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine ranlux_seeds(flag,index,rvalue)
c     ******************************************************************
c     * the seeds global variable in the random number generator       *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag,index
      real*8 rvalue
c     --- information variables
      real*8 seeds
      dimension seeds(24)
      save seeds

      if ( flag .eq. GLB_SET ) then
         seeds(index) = rvalue
      elseif ( flag .eq. GLB_GET ) then
         rvalue = seeds(index)
      else
         write(6,*) 'RANLUX_seeds: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine ranlux_twomtwelve(flag,rvalue)
c     ******************************************************************
c     * the twomtwelve global variable in the random number generator  *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag
      real*8 rvalue
c     --- information variables
      real*8 twomtwelve
      save twomtwelve

      if ( flag .eq. GLB_SET ) then
         twomtwelve = rvalue
      elseif ( flag .eq. GLB_GET ) then
         rvalue = twomtwelve
      else
         write(6,*) 'RANLUX_twomtwelve: invalid global flag',flag
         stop
      endif

      return
      end
      
      subroutine ranlux_twomtwentyfour(flag,rvalue)
c     ******************************************************************
c     * the twomtwentyfour global variable in the random number        *
c     * generator                                                      *
c     *                                                                *
c     * originally written 11-10-2006 by M.G. Martin                   *
c     * last modified 11-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer*4 flag
      real*8 rvalue
c     --- information variables
      real*8 twomtwentyfour
      save twomtwentyfour

      if ( flag .eq. GLB_SET ) then
         twomtwentyfour = rvalue
      elseif ( flag .eq. GLB_GET ) then
         rvalue = twomtwentyfour
      else
         write(6,*) 'RANLUX_twomtwentyfour: invalid global flag',flag
         stop
      endif

      return
      end

      function rngdebug()
c     ******************************************************************
c     * a not at all random option for the RNG used just for debugging *
c     *                                                                *
c     * originally written 02-05-2013 by M.G. Martin                   *
c     * last modified 02-05-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
      double precision rngdebug

      rngdebug = 0.3d0

      return
      end
