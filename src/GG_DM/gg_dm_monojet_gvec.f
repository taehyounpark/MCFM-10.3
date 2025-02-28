!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_dm_monojet_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'qcdcouple.f'
c  in is the label of the momentum contracted with n
      integer:: j,k,in
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: n(4),p(mxpart,4),hdecay,fac,
     & qqgDMn,gggDMn,p1p2(-1:1,-1:1)

      msq(:,:)=0._dp


      call dmsdecay(p,3,4,hdecay)
c      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      fac=hdecay/dm_lam**6*as**2*16d0

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

      if (in == 1) then
      p1p2(0,-1)=-aveqg*fac*qqgDMn(2,5,1,p,n)
      p1p2(0,+1)=-aveqg*fac*qqgDMn(2,5,1,p,n)
      p1p2(0,0)=+avegg*fac*gggDMn(5,2,1,p,n)
      elseif (in == 2) then
      p1p2(+1,0)=-aveqg*fac*qqgDMn(1,5,2,p,n)
      p1p2(-1,0)=-aveqg*fac*qqgDMn(5,1,2,p,n)
      p1p2(0,0)=+avegg*fac*gggDMn(1,5,2,p,n)
      elseif (in == 5) then
      p1p2(1,-1)=+aveqq*fac*qqgDMn(1,2,5,p,n)
      p1p2(-1,1)=+aveqq*fac*qqgDMn(2,1,5,p,n)
      p1p2(0,0)=+avegg*fac*gggDMn(1,2,5,p,n)
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k == -j)) then
          msq(j,k)=p1p2(1,-1)
      elseif ((j < 0) .and. (k == -j)) then
          msq(j,k)=p1p2(-1,1)
      elseif ((j == 0) .and. (k == 0)) then
          msq(j,k)=p1p2(0,0)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &    p1p2(+1,0)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    p1p2(-1,0)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    p1p2(0,+1)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    p1p2(0,-1)
      endif
      enddo
      enddo

      return
      end

      function qqgDMn(j1,j2,j5,p,n)
      implicit none
      include 'types.f'
      real(dp):: qqgDMn

c---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> H((p3+p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'

      integer:: j1,j2,j5
      real(dp):: Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,dot,s,t,u

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)


c     Asq=(as/(3d0*pi))**2/vevsq
c      Asq=one/dm_lam**6
      Asq=one

      s=2d0*Dot(p,j1,j2)
      t=2d0*Dot(p,j1,j5)
      u=2d0*Dot(p,j2,j5)

c--- RKE Schoonship answer, ncalc.out
c     Id,bit=g^2*A^2*V/2*(2*(nDp1*u-nDp2*t)^2/s^2
c           + 0.5*nDn*(u+t)^2/s)

      qqgDMn=-Asq*gsq*V/2d0*(2d0*(nDp1*u-nDp2*t)**2/s**2
     &                    +0.5d0*nDn*(u+t)**2/s)

      return
      end

      function gggDMn(j1,j2,j5,p,n)
      implicit none
      include 'types.f'
      real(dp):: gggDMn

c---calculates the amplitude squared for the process
c   g(p1)+g(p2) --> H((p3+p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'

      integer:: j1,j2,j5
      real(dp):: Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,dot,s,t,u,sh

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)

c      Asq=(as/(3d0*pi))**2/vevsq
c     Asq=one/dm_lam**6
      Asq=one

      s=2d0*Dot(p,j1,j2)
      t=2d0*Dot(p,j1,j5)
      u=2d0*Dot(p,j2,j5)
      sh=s+t+u

c--- JMC answer, gggH.frm
c -f(b,a,c)^2/s12/s13/s23*(
c -(n.n)/2*(s12^4+s13^4+s23^4+mHsq^4-2*(s13^2*s23^2+s12^2*mHsq^2))
c +2*(p1.n*s23-p2.n*s13)^2/s13/s23*(s13^2*s23^2+s12^2*mHsq^2)/s12);

      gggDMn=Asq*gsq*V*xn*(
     & -nDn/2d0*(s**4+t**4+u**4+sh**4-2d0*(t**2*u**2+s**2*sh**2))
     & +2d0*(nDp1*u-nDp2*t)**2/t/u*(t**2*u**2+s**2*sh**2)/s)/s/t/u

      return
      end
