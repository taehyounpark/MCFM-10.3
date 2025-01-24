!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function virtgamgam_msbar(s,t,u)
      implicit none
      include 'types.f'
      real(dp):: virtgamgam_msbar

      !include 'epinv.f'
      !include 'epinv2.f'
      include 'constants.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'scet_const.f'
      real(dp):: s,t,u
      complex(dp):: m1(2,2,2,2),m2fin(2,2,2,2),lnrat,
     & I1ggtogamgam,xlog,FL(2,2,2,2),FSL(2,2,2,2)
      complex(dp) :: I1toZ
      real(dp), parameter :: epinv = 0._dp, epinv2 = 0._dp
      include 'scet_beta.f'

      xlog=lnrat(musq,-s)
c--- This is taken from hep-ph/0109078 Eq.(2.11); note however that the log
c--- proportional to the beta-function coefficient is added in Eq. (2.11)
c--- and subtracted again in Eq. (4.5), therefore we omit it here.
      I1ggtogamgam=
     & -xn*((epinv*epinv2+epinv*xlog+0.5_dp*xlog**2) - pisq/12._dp
     & +(11._dp/6._dp-real(nflav,dp)/(3._dp*xn))*epinv)

      I1toZ = -xn/2._dp*(xlog**2 + be0/xn*xlog) + pisq/12._dp*xn

      call M1fill(s,t,u,m1)
      call twoloop(s,t,u,FL,FSL)

      m2fin = xn*FL - FSL/xn
      virtgamgam_msbar = sum(2._dp*ason2pi
     & *real(conjg(m1)*(m1*I1ggtogamgam -m1*I1toZ + m2fin)))

      end function
