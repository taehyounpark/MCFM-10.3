!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine phi1_2m_nobw(m2,x3,xth,xphi,s3min,p1,p2,p3,wt,*)
      implicit none
      include 'types.f'
c     massive particle p1 decaying into p2 mass m2 and p3 mass-squared s3.
c     with invariant mass of particle three s3 integrated over.
c     vectors returned p2 and p3 are in the same frame as p1 is supplied
c     Expression evaluate is
c     ds2 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6
c     delta(p2^2-s2) delta(p3^2-s3)
      include 'constants.f'
      include 'breit.f'
      real(dp):: p1(4),p2(4),p3(4),p3cm(4)
      real(dp):: x3,xth,xphi,costh,sinth,phi,cphi,sphi
      real(dp):: wt,wt0,w3
      real(dp):: s3max,s3min
      real(dp):: m1,m2,m3,s1,s2,s3,lambda
      integer:: j
      parameter(wt0=one/8._dp/pi)

      wt=0._dp
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      if (s1 <= 0._dp) return 1
      m1=sqrt(s1)
      s2=m2**2
      s3max=(m2-m1)**2
      if (s3min > s3max) return 1
      w3=s3max-s3min
      s3=s3max*x3+s3min*(1._dp-x3)

      if (s3 <= 0._dp) return 1
      m3=sqrt(s3)
      costh=two*xth-one
      phi=twopi*xphi
      sinth=sqrt(max(one-costh**2,zip))
      cphi=cos(phi)
      sphi=sin(phi)
      lambda=((s1-s2-s3)**2-4._dp*s2*s3)
      if (lambda < 0._dp) then
      write(6,*) 'lambda in phi1_2m',lambda
      write(6,*) 's1 in phi1_2m',s1
      write(6,*) 's2 in phi1_2m',s2
      write(6,*) 's3 in phi1_2m',s3
      write(6,*) 'm1 in phi1_2m',m1
      write(6,*) 'm2 in phi1_2m',m2
      write(6,*) 'm3 in phi1_2m',m3
      write(6,*) 'x3 in phi1_2m',x3
      write(6,*) 'mass3 in phi1_2m',mass3
      return 1
      endif
      lambda=sqrt(lambda)
      wt=wt0*w3*lambda/s1

      p3cm(4)=m1/two*(s1+s3-s2)/s1
      p3cm(1)=m1/two*lambda/s1*sinth*sphi
      p3cm(2)=m1/two*lambda/s1*sinth*cphi
      p3cm(3)=m1/two*lambda/s1*costh
      call boost(m1,p1,p3cm,p3)
      do j=1,4
      p2(j)=p1(j)-p3(j)
      enddo
      if (  (p1(4) < 0._dp)
     & .or. (p2(4) < 0._dp)
     & .or. (p3(4) < 0._dp)) then
c      write(6,*) 'p1(4)',p1(4)
c      write(6,*) 'p2(4)',p2(4)
c      write(6,*) 'p3(4)',p3(4)
c      write(6,*) 'p1sq',p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2,s1
c      write(6,*) 'p2sq',p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2,s2
c      write(6,*) 'p3sq',p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2,s3
c      write(6,*) 'in phi1_2m.f'
c      write(6,*) n2,n3
        return 1
      endif
      return
      end



