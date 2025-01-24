!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function Hqarbsq(i1,i2,i3,i4)
      implicit none
      include 'types.f'
      real(dp):: Hqarbsq

c----Matrix element squared for the process
c----q(p1)+r(p3) -> q(p2)+r(p4)+Higgs
c----summed over incoming/outgoing colors and spins.
c    with a factor of gsq*Asq removed.

c       +V*((s13*s24-s23*s14)**2+s12**2*s34**2)/s34**2/s12**2
c       +1/2*V*((s13-s24)**2+(s14-s23)**2)/s34/s12
c       -1/2*e*V*(s13+s14+s23+s24)**2/s34/s12;

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
c GZ
c      write(*,*) 'in Hqarbsq'
c       s12 = 1.000000000000000_dp
c       s13= -0.2861267766967376_dp
c       s14= -0.3784302994364646_dp
c       s23= -0.6215181025806800_dp
c       s24= -0.3421507931851938_dp
c       s34=  0.6319744100498450_dp

      Hqarbsq=V*(((s13*s24-s23*s14)**2+s12**2*s34**2)/s34**2/s12**2
     & +0.5_dp*((s13-s24)**2+(s14-s23)**2)/s34/s12)
c      write(*,*) 'In Hqarbsq',Hqarbsq
      return
      end
