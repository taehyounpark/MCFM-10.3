!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Hqaqasq(i1,i2,i3,i4)
      implicit none
      include 'types.f'
      real(dp):: Hqaqasq

c----Interference piece of Matrix element squared for the process
c----q(p1)+q(p3) -> q(p2)+q(p4)+Higgs
c----with identical quarks
c----summed over incoming/outgoing colors and spins.
c    with a factor of gsq*Asq removed.

c      Cf*((s13-s24)^2*(s12*s34-s13*s24+s14*s23)
c         -2*(-s12*s34+s13*s24+s14*s23)*(s12*s34+s13*s24-s14*s23))
c          /(s12*s14*s23*s34)
c       +Cf*e*(s12*s34-s13*s24+s14*s23)
c        *((s23+s14)*(s12+s34)-(s13+s24)^2)
c         /(s12*s14*s23*s34)
c       -Cf*e^2*(s24+s23+s14+s13)*(s34+s24+s13+s12)*(s12*s34-s13*s24+s14*s23)
c       /(s12*s14*s23*s34);

      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4
      real(dp):: s12,s13,s14,s23,s24,s34
      s12=s(i1,i2)
      s13=s(i1,i3)
      s14=s(i1,i4)
      s23=s(i2,i3)
      s24=s(i2,i4)
      s34=s(i3,i4)

      Hqaqasq=+Cf*((s13-s24)**2*(s12*s34-s13*s24+s14*s23)
     &    -2._dp*(-s12*s34+s13*s24+s14*s23)*(s12*s34+s13*s24-s14*s23))
     &     /(s12*s14*s23*s34)

      return
      end
