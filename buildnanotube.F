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
      subroutine twh_buildnanotube(lfailure,imolty)
c     ******************************************************************
c     * this subroutine was submitted to the Towhee project by Steve   *
c     * Lustig and he had adapted it from a previous routine as        *
c     * described here with the original comments                      *
c     * --- start of the original comments --                          *
c     * adaptations made by Steve Lustig                               *
c     * 0) corrects several nonANSI aspects for successful local       *
c     *    compilation                                                 *
c     * 1) combines a-2 and a-3 to provide extended unit cells         *
c     * 2) writes pdb file of results                                  *
c     * version 2001.11.15                                             *
c     * Atomic Coordinates in the unit cell of a (n,m) nanotube        *
c     * Made by 10/31/95   by T. Takeya and R. Saito                   *
c     * Reference: ``Physical Properties of Carbon Nanotubes''         *
c     *            Imperial College Press, by R. Saito, G. Dresselhaus,*
c     *            and M.S. Dresselhaus                                *
c     * Input: n,m  (from keyboard)                                    *
c     * Output: tube.xyz  (coordinates. See detail in the above book.) *
c     * en.xyz2   (a work file for the program A-3)                    *
c     * --- end of the original comments  ---                          *
c     * many things have now changed and the original comments are no  *
c     * longer valid.                                                  *
c     *                                                                *
c     * routine originally incorporated on 05-17-2004 by M.G. Martin   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty
c     --- local variables
      logical lautotorsion
c     --- character string scalars
      character*10 forcefield
      character*15 idname
      character*30 filename
c     --- character string arrays
      character*10 atomname
      dimension atomname(NUMAX)
c     --- integer scalars
      integer n,m,ncells,totatoms
      integer np,nq,ndr,nn,ii,i,k,kk2,kk,j,ivib,iunit,junit
c     --- double precision scalars
      double precision dz,bondlength,qtemp
      double precision sq3,pi,a,r,c,t,rs,q1,q2,q3,q4,q5,h1,h2,x1,y1,z1
      double precision z3,x2,y2,z2,test,rx,ry,rz,rij
c     --- double precision arrays
      double precision x,y,z
      dimension x(NUMAX),y(NUMAX),z(NUMAX) 

      idname = 'buildnanotube'

      write(6,'(a51,i5)') 
     &     'Building the nanotube input file for molecule type:',imolty
c     --- read information from towhee_input
c     --- forcefield
      call twh_read_labeled_string(lfailure,4,'forcefield',.true.,idname
     &     ,forcefield)
      if ( lfailure ) return
c     --- read in the repeat atomname
      call twh_read_labeled_string(lfailure,4,'atomname',.true.,idname
     &     ,atomname(1))
      if ( lfailure ) return
c     --- read in the qqatom
      qtemp = twh_read_labeled_float(lfailure,4,'qqatom',.false.,
     &    idname)
      if ( lfailure ) return
      call twh_qqatom(GLB_SET,imolty,1,qtemp)
      write(6,'(a7,f10.5)') 'qqatom:',qtemp

c     --- read in the 'n' parameter
      n = twh_read_labeled_integer(lfailure,4,'nanotube_n',.true.,
     &    idname)
      if ( lfailure ) return
      if ( n .lt. 0 ) then
         write(6,*) 'BUILDNANOTUBE: nanotube_n must be non-negative'
         lfailure = .true.
         return
      endif

c     --- read in the 'm' parameter
      m = twh_read_labeled_integer(lfailure,4,'nanotube_m',.true.,
     &    idname)
      if ( lfailure ) return
      if ( m .lt. 0 ) then
         write(6,*) 'BUILDNANOTUBE: nanotube_m must be non-negative'
         lfailure = .true.
         return
      endif

c     --- read in the number of repeat cells
      ncells = twh_read_labeled_integer(lfailure,4,'nanotube_ncells'
     &     ,.true.,idname)
      if ( lfailure ) return
      if ( ncells .lt. 1 ) then 
         write(6,*) 'BUILDNANOTUBE: nanotube_ncells must be positive'
         lfailure = .true.
         return
      endif

c     --- read in the atom bond length
      bondlength = twh_read_labeled_float(lfailure,4
     &     ,'nanotube_bondlength',.false.,idname)
      if ( lfailure ) return
      write(6,'(a20,f10.5)') 'nanotube_bondlength:',bondlength
      if ( bondlength .le. 0.0d0 ) then
         write(6,*) 'BUILDNANOTUBE: nanotube_bondlength must be postive'
         lfailure = .true.
         return
      endif

c     --- eventually will put error check stuff here once I figure 
c     --- out the proper values for n and m MGM

      call gen11(lfailure,n,m,np,nq,ndr) 
      if ( lfailure ) return
c     a: the length of the unit vector 
c     pi = 3.141592... sq3 = 1.732.... 
      sq3=sqrt(3.0d0) 
      pi=4.0d0*atan(1.0d0) 
      a=sqrt(3.0d0)*bondlength 
c     r;|R| , c;|C_h| , t;|T| 
      r=a*sqrt(dfloat(np*np+nq*nq+np*nq)) 
      c=a*sqrt(dfloat(n*n+m*m+n*m)) 
      t=sqrt(3.0d0)*c/dfloat(ndr) 
      write(6,'(a26,f10.5)') 'nanotube unit cell length:',t 
      write(6,'(a16,f10.5)') 'nanotube radius:',c/2.0d0/pi 
c     nn: the number of hexagon in the unit cell N 
      nn=2*(n**2+m**2+m*n)/ndr 
      write(6,'(a25,i5)') 'nanotube unit cell atoms:',nn*2 
c     rs: radius of the tube 
      if (2*nn .gt. NUMAX ) then
         write(6,*) 'parameter NUMAX is too small' 
         lfailure = .true.
         return
      endif
      rs=c/(2.0d0*pi) 
c     q1: the chiral angle for C_h 
      q1=atan((sq3*dfloat(m))/dfloat(2*n+m)) 
c     q2: the 'chiral' angle for R 
      q2=atan((sq3*dfloat(nq))/dfloat(2*np+nq)) 
c     q3: the angle between C_h and R 
      q3=q1-q2 
c     q4: a period of an angle for the A atom 
      q4=2.0d0*pi/dfloat(nn) 
c     q5: the difference of the angle between the A and B atoms 
      q5=bondlength*cos((pi/6.0d0)-q1)/c*2.0d0*pi 
c     h1: 
      h1=abs(t)/abs(sin(q3)) 
c     h2: Delta z between the A and B atoms 
      h2=bondlength*sin((pi/6.0d0)-q1) 
C     write(6,*) 'q1: =',q1*180.0d0/pi 
C     write(6,*) 'q2: =',q2*180.0d0/pi 
C     write(6,*) 'q4: =',q4*180.0d0/pi 
C     write(6,*) 'q5: =',q5*180.0d0/pi 
c     Calculate 2*nn atoms in the unit cell. 
c     The A atom 
      ii=0 
      do i=0,nn-1 
         x1=0 
         y1=0 
         z1=0 
         k=int(dfloat(i)*abs(r)/h1) 
         x1=rs*cos(dfloat(i)*q4) 
         y1=rs*sin(dfloat(i)*q4) 
         z1=(dfloat(i)*abs(r)-dfloat(k)*h1)*sin(q3) 
         kk2=abs(int(z1/t))+1 
c        Check the A atom is in the unit cell 0 < z1 < t 
         if(z1.gt.t-0.02)then 
            z1=z1-t*dfloat(kk2) 
         endif 
         if(z1.lt.-0.02) then 
            z1=z1+t*dfloat(kk2) 
         endif 
         ii=ii+1 
         x(ii)=x1 
         y(ii)=y1 
         z(ii)=z1 

c        The B atoms 
         z3=(dfloat(i)*abs(r)-dfloat(k)*h1)*sin(q3)-h2 
         ii=ii+1 
c        Check the B atom is in the unit cell 0 < z3 < t 
         if((z3.ge.-0.02).and.(z3.le.t-0.02))then 
c           yes 
            x2 =rs*cos(dfloat(i)*q4+q5) 
            y2 =rs*sin(dfloat(i)*q4+q5) 
c            z2 =dfloat(dfloat(i)*abs(r)-dfloat(k)*h1)*sin(q3)-h2 
            z2 =(dfloat(i)*abs(r)-dfloat(k)*h1)*sin(q3)-h2 
            x(ii)=x2 
            y(ii)=y2 
            z(ii)=z2 
         else 
c           no 
            x2 =rs*cos(dfloat(i)*q4+q5) 
            y2 =rs*sin(dfloat(i)*q4+q5) 
            z2 =(dfloat(i)*abs(r)-dfloat(k+1)*h1)*sin(q3)-h2 
            kk=abs(int(z2/t))+1 
            if(z2.gt.t-0.01)then 
               z2=z2-t*dfloat(kk) 
            endif 
            if(z2.lt.-0.01) then 
               z2=z2+t*dfloat(kk) 
            endif 
            x(ii)=x2 
            y(ii)=y2 
            z(ii)=z2 
         endif 
      enddo

c     --- compute the true number of total atoms and assign their
c     --- positions
      totatoms = 0
      do j=0,ncells-1 
         dz=j*t 
         do i=1,nn*2 
            totatoms = totatoms + 1
c           --- safety check totatoms
            if ( totatoms .gt. NUMAX ) then
               write(6,*) 'BUILDNANOTUBE: totatoms exceeds NUMAX'
               lfailure = .true.
               return
            endif
            x(totatoms) = x(i)
            y(totatoms) = y(i)
            z(totatoms) = z(i) + dz
         enddo 
      enddo 
      write(6,'(a21,i5)') 'nanotube total atoms:',totatoms
c     Output to the file towhee_nanotube
      filename = 'towhee_nanotube'
      call twh_openfile(lfailure,60,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return
      write(60,*) 'nanotube structure',totatoms
      do i=1,totatoms
         write(60,'(3f10.5)')x(i),y(i),z(i) 
      enddo 
      close(60) 

c     --- setup logic for assemble
c     --- assign the default atomname and qqatom to all units
      call twh_qqatom(GLB_GET,imolty,1,qtemp)
      do iunit = 2,totatoms
         atomname(iunit) = atomname(1)
         call twh_qqatom(GLB_SET,imolty,iunit,qtemp)
      enddo
c     --- loop through all atoms and detect bond neighbors 
c     --- as those within the specified bondlength
c     --- no need for the minimum image convention as this is a single
c     --- unit cell template
      test = 1.01d0*bondlength
      do iunit = 1,totatoms-1
         do junit = iunit+1,totatoms
            rx = x(junit)-x(iunit)
            ry = y(junit)-y(iunit)
            rz = z(junit)-z(iunit)
            rij = twh_distance(rx,ry,rz)
            if ( rij .lt. test ) then
c              --- bond neighbor found, add to both data structures
               call twh_invib(GLB_INCR,imolty,iunit,1)
               call twh_invib(GLB_GET,imolty,iunit,ivib)
               call twh_ijvib(GLB_SET,imolty,iunit,ivib,junit)
               call twh_invib(GLB_INCR,imolty,junit,1)
               call twh_invib(GLB_GET,imolty,junit,ivib)
               call twh_ijvib(GLB_SET,imolty,junit,ivib,iunit)
            endif
         enddo
      enddo
c     --- autodetect all torsions
      lautotorsion = .true.
      call twh_assemble(lfailure,imolty,totatoms,0,forcefield,atomname
     &     ,'standard','manual',lautotorsion)
      if ( lfailure ) return

c     --- set nmaxcbmc
      call twh_nmaxcbmc(GLB_SET,imolty,totatoms)

      return
      end 

      subroutine gen11(lfailure,n,m,np,nq,ndr) 
c     ******************************************************************
c     * original comment: 10/16/95   R. Saito                          *
c     *                                                                *
c     * incorporated into Towhee 05-17-2004 by M.G. Martin             *
c     * last modified 04-21-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- functions local to buildnanotube
      logical lfailure
      integer igcm
c     --- variables passed to/from the subroutine
      integer n,m,np,nq,ndr
c     --- local variables
      integer nnp,nnq
      dimension nnp(100),nnq(100) 
      integer itest,itest1,nd,l2,l,dt,nr,ns,nt2,nt,nn,ichk,n60,itest2
      integer j2,j1,itest3,i

      double precision a,eps

      itest=1 
      itest1=0 
      nd=igcm(n,m) 
      if(mod((n-m),3*nd).eq.0) then 
         ndr=3*nd 
      else 
         ndr=nd 
      endif 
      if(itest.eq.1) nd=nd 
      if(itest.eq.1) ndr=ndr 
      a=sqrt(3.0)*1.421 
      eps=1.0e-5 
c     L 
      l2=n*n+m*m+n*m 
      if(l2.le.0) then
         write(6,*) 'GEN11: error l2.le.0'
         lfailure = .true.
         return
      endif
      l = int(sqrt(dfloat(l2))+eps) 
      dt = int(a*sqrt(dfloat(l2))/3.1415926525)
c     T 
      nr=(2*m+n)/ndr 
      ns=-(2*n+m)/ndr 
      nt2=3*l2/ndr/ndr 
      nt=int(sqrt(dfloat(nt2))+eps) 
c     N 
      nn=2*l2/ndr 
c     R 
      ichk=0 
      if(nr.eq.0) then 
         n60=1 
      else 
         n60=nr 
      endif 
      itest2=1 
      do np=-abs(n60),abs(n60) 
         do nq=-abs(ns),abs(ns) 
            j2 = nr*nq - ns*np 
            if(j2.eq.1) then 
               j1 = m *np -  n*nq 
               if( j1.gt.0 .and. j1.lt.nn ) then 
                  ichk=ichk+1 
                  nnp(ichk)=np 
                  nnq(ichk)=nq 
               endif 
            endif 
         enddo
      enddo
      if(ichk.eq.0) then 
         write(6,*) 'GEN11: not found p,q strange!!' 
         lfailure = .true.
         return
      endif 
      itest3=1 
      if(ichk.ge.2) then 
         if(itest3.eq.1) then 
            write(6,*) 'GEN11: more than 1 pair of p,q strange!!' 
            lfailure = .true.
            return
         endif 
         if(nr.ne.0 .and. ns.ne.0) then 
            if(itest1.eq.1) then 
               do i=1,ichk 
                  if((m*nnp(i)-n*nnq(i)).lt.nn) goto 777 
               enddo
            endif 
         endif 
      endif 
 777  continue 
      if(itest.eq.1) then 
         np=nnp(1) 
         nq=nnq(1) 
      endif 
      return 
      lfailure = .true.
      return
      end


      function igcm(ii,jj) 
c     ******************************************************************
c     * calculate the highest common divisor                           *
c     * originally added into Towhee 05-17-2004 by M.G. Martin         *
c     * last modified 06-17-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- functions local to buildnanotube
      integer igcm
c     --- variables passed to/from the function
      integer ii,jj
c     --- local variables
      integer i,j,iw,ir

      i=abs(ii) 
      j=abs(jj) 
      if(j.gt.i) then 
         iw=j 
         j=i 
         i=iw 
      endif 
      if(j.eq.0) then 
         igcm=i 
         return 
      endif 
 10   continue 
      ir=mod(i,j) 
      if(ir.eq.0) then 
         igcm=j 
         return 
      else 
         i=j 
         j=ir 
      endif 
      goto 10 
      return
      end 
