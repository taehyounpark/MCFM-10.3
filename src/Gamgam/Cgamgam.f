!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function Cgamgam(i1,i2,i3,i4,i5,i)
      implicit none
      include 'types.f'
      real(dp):: Cgamgam

cCCCCC Matrix element squared for
c     qbar(-p1)+q(-p2)=gamma(p3)+gamma(p4)+g(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'zcouple_cms.f'
      integer:: i1,i2,i3,i4,i5,i
      real(dp),parameter::statfac=0.5_dp
      Cgamgam=s(i1,i2)
     &*((s(i2,i5)**2+s(i1,i5)**2)/(s(i1,i4)*s(i2,i4)*s(i1,i3)*s(i2,i3))
     & +(s(i2,i4)**2+s(i1,i4)**2)/(s(i1,i5)*s(i2,i5)*s(i1,i3)*s(i2,i3))
     & +(s(i2,i3)**2+s(i1,i3)**2)/(s(i1,i5)*s(i2,i5)*s(i1,i4)*s(i2,i4)))

      Cgamgam=16._dp*cf*xn*gsq*abs(zesq)**2*Q(i)**4*Cgamgam*statfac

      return
      end
