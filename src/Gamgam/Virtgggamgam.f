!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function virtgamgam(s,t,u)
      implicit none
      include 'types.f'
      real(dp):: virtgamgam

      include 'epinv.f'
      include 'epinv2.f'
      include 'constants.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'nflav.f'
      real(dp):: s,t,u
      complex(dp):: m1(2,2,2,2),m2fin(2,2,2,2),lnrat,
     & I1ggtogamgam,xlog,FL(2,2,2,2),FSL(2,2,2,2)

      xlog=lnrat(musq,-s)
c--- This is taken from hep-ph/0109078 Eq.(2.11); note however that the log
c--- proportional to the beta-function coefficient is added in Eq. (2.11)
c--- and subtracted again in Eq. (4.5), therefore we omit it here.
      I1ggtogamgam=
     & -xn*((epinv*epinv2+epinv*xlog+0.5_dp*xlog**2)
     & +(11._dp/6._dp-real(nflav,dp)/(3._dp*xn))*epinv)

      call M1fill(s,t,u,m1)

c--- testing the 1-loop matrix element contribution up to Order(ep^2)
c      call oneloopep(s,t,u,IxM1ep)

c--- pure 2-loop contribution
      call twoloop(s,t,u,FL,FSL)

c--- Compare m1 with m1ep
c      do h1=1,2
c      do h2=1,2
c      do h3=1,2
c      do h4=1,2
c      write(6,*) h1,h2,h3,h4,m1(h1,h2,h3,h4)*I1ggtogamgam,
c     &           IxM1ep(h1,h2,h3,h4)
c      enddo
c      enddo
c      enddo
c      enddo
c      pause

      m2fin = xn*FL - FSL/xn
      virtgamgam = sum(2._dp*ason2pi*real(conjg(m1)*(m1*I1ggtogamgam + m2fin)))

      end
