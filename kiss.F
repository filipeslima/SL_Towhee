c     Return double precision pseudo-random number in [0,1) on each call.
c
      function kiss99()
      implicit none
      integer*4 result
      real*8 kiss99
      call kiss(result)
c     the constants in the following statement are 2**31 and 2**32
      kiss99=(result+2147483648D0)/4294967296D0
      end


c     Marsaglia's KISS (Keep It Simple Stupid) random number generator.
c     Combines:
c        (1) The congruential generator x(n)=69069*x(n-1)+1327217885,
c        period 2^32,
c        (2) A 3-shift shift-register generator, period 2^32-1,
c        (3) Two 16-bit multiply-with-carry generators, period
c        597273182964842497>2^59.  Overall period>2^123.
c     Default seeds x,y,z,w.
c     Get random integer with kiss(); set your own seeds with kisset.
c
      subroutine kiss(result)
      implicit none
      integer*4 result, x,y,z,w, iseed, m,k,n, nrng,irng(*)
      save x,y,z,w
      data x,y,z,w/123456789,362436069,521288629,916191069/
c     Statement function
      m(k,n) = ieor (k, ishft (k, n) )
c
      x = 69069 * x + 1327217885
      y = m (m (m (y, 13), - 17), 5)
      z = 18000 * iand (z, 65535) + ishft (z, - 16)
      w = 30903 * iand (w, 65535) + ishft (w, - 16)
      result = x + y + ishft (z, 16) + w
      return

C     Entry point for kiss initialization
      entry kissinit(iseed)
      x = iseed
      y = iseed
      z = iseed
      w = iseed
      return

C     Entry point for kiss snapshot
      entry kisssave(nrng,irng)
      nrng = 4
      irng(1) = x
      irng(2) = y
      irng(3) = z
      irng(4) = w
      return

C     Entry point for kiss restart (re-initialization)
      entry kissrestart(nrng,irng)
      if (nrng.ne.4) then
         write(6,*) 'In kissrestart, nrng has wrong value: ', nrng
         stop
      endif
      x = irng(1)
      y = irng(2)
      z = irng(3)
      w = irng(4)
      return

C     All access subroutines defined
C     Complete subroutine definition
      end

