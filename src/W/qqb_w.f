!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_w(p,msq)
      implicit none
      include 'types.f'

c----Matrix element for W production
c----averaged over initial colours and spins
c For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4))
c For nwz=-1
c     d(-p1)+ubar(-p2)-->W^-(e^-(p3)+nbar(p4))
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zcouple_cms.f'
      include 'ckm.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qqb,qbq,s

c--statement function
      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))


      fac=abs((zesq/zxw)**2)*xn
c--   calculate propagator
      fac=aveqq*fac/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
c---case dbar-u or ubar-d
      qqb=fac*s(1,4)**2
      qbq=fac*s(2,4)**2

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=zip
          if ((j > 0) .and. (k < 0)) then
            msq(j,k)=Vsq(j,k)*qqb
          elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=Vsq(j,k)*qbq
          endif
      enddo
      enddo
      return
      end
