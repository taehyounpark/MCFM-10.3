!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine w2jetsq_mass(i1,i2,i3,i4,i5,i6,p,msq)
      implicit none
      include 'types.f'

c     s(-p1)+cbar(-p2) --> l(p3)+abar(p4)+g(p5)+g(p6)
c     with c massive
c     multiplied by (((a+l)^2-M**2)/(a+l)^2)^2*g^4/gwsq^2/2
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'mmsq_cs.f'
      complex(dp):: qcd1(2,2,2),qcd2(2,2,2),qed(2,2,2)
      real(dp):: msq1,msq2,msqq,msq,p(mxpart,4),q(mxpart,4),a5,dot
      real(dp):: uc,us,dc,ds,cs,p156,p256,mch,
     & invtwog1Dc,invtwog2Dc
      integer:: i1,i2,i3,i4,i5,i6,nu,j,h1,h2,hf

      uc=dot(p,i1,i5)
      us=dot(p,i1,i6)
      dc=dot(p,i2,i5)
      ds=dot(p,i2,i6)
      cs=dot(p,i5,i6)

      mch=sqrt(abs(p(i2,4)**2-p(i2,1)**2-p(i2,2)**2-p(i2,3)**2))
      p156=2._dp*(uc+us+cs)
      p256=2._dp*(dc+ds+cs)

c----define massless momenta
      a5=0.5_dp*mch**2/dot(p,i1,i2)

      do nu=1,4
      do j=1,6
        if (j==i2) then
          q(j,nu)=p(i2,nu)-p(i1,nu)*a5
        else
          q(j,nu)=p(j,nu)
        endif
      enddo
      enddo

      invtwog1Dc=1._dp/2._dp/dot(p,i5,i2)
      invtwog2Dc=1._dp/2._dp/dot(p,i6,i2)

      call spinoru(6,q,za,zb)

      call subqcdm(i1,i2,i3,i4,i5,i6,p156,p256,za,zb,
     & invtwog1Dc,invtwog2Dc,mch,qcd1,qcd2)

      do hf=1,2
      qed(2,2,hf)=qcd1(2,2,hf)+qcd2(2,2,hf)
      qed(2,1,hf)=qcd1(2,1,hf)+qcd2(2,1,hf)
      qed(1,2,hf)=qcd1(1,2,hf)+qcd2(1,2,hf)
      qed(1,1,hf)=qcd1(1,1,hf)+qcd2(1,1,hf)
      enddo

      msq1=0._dp
      msq2=0._dp
      msqq=0._dp
      do h1=1,2
      do h2=1,2
      do hf=1,2
      msq1=msq1+abs(qcd1(h1,h2,hf))**2
      msq2=msq2+abs(qcd2(h1,h2,hf))**2
      msqq=msqq+abs(qed(h1,h2,hf))**2
      enddo
      enddo
      enddo

      mmsq_cs(0,1,1)=-ninth*msqq
      mmsq_cs(1,1,1)=msq1
      mmsq_cs(2,1,1)=msq2

      msq=mmsq_cs(0,1,1)+mmsq_cs(1,1,1)+mmsq_cs(2,1,1)

      return
      end
