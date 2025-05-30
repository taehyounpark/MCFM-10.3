!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      !complex(dp) w(35)

      w(1)=s**(-1)
      w(2)=Log(mt2*s**(-1))
      w(3)=Hr1(0)
      w(4)=Hr1(1)
      w(5)=Hr2(0,0)
      w(6)=Hr2(0,1)
      w(7)=Hr2(1,0)
      w(8)=Hr2(1,1)
      w(9)=1/( - 1 + u2)
      w(10)=1/( - 2*u2 + 2*u2**2)
      w(11)=1/( - u2 + u2**2)
      w(12)=1/( - u2 + 3*u2**2 - 3*u2**3 + u2**4)
      w(13)=1/(u2 - 2*u2**2 + u2**3)
      w(14)=u2**(-1)
      w(15)=1/(s - 2*s*u2 + s*u2**2)
      w(16)=1/(s*u2 - 2*s*u2**2 + s*u2**3)
      w(17)=1/(s*u2**2 - 2*s*u2**3 + s*u2**4)
      w(18)=w(6) + w(7)
      w(19)=Pi**2
      w(20)=im*Pi
      w(21)=w(20) + w(2)
      w(22)=w(3)*w(21)
      w(22)= - w(18) + w(22) + w(19) + w(5)
      w(23)=2*w(8)
      w(24)=2*w(20)
      w(25)= - w(24) - 3*w(2)
      w(25)=w(2)*w(25)
      w(22)=w(23) + w(25) + 2*w(22)
      w(22)=w(1)*w(22)
      w(21)=w(21)*w(4)
      w(25)=2*w(1)
      w(25)=w(25)*w(21)
      w(22)=w(25) + w(22)
      w(22)=2*w(22)
      w(25)=4*w(20)
      w(26)= - w(25) + w(19) - 12
      w(27)=u2*w(26)
      w(28)=8*w(20)
      w(27)=w(27) + w(28) + 24 - w(19)
      w(27)=u2*w(27)
      w(27)= - 24 + w(27)
      w(27)=u2*w(27)
      w(27)=12 + w(27)
      w(27)=w(13)*w(27)
      w(21)=w(21) + w(8) - w(18)
      w(21)=w(9)*w(21)
      w(29)=2*u2
      w(30)= - 2 - w(20)
      w(30)=w(30)*w(29)
      w(30)=w(30) + 4 + w(20)
      w(30)=u2*w(30)
      w(31)=u2 - 2
      w(32)= - w(3)*w(31)
      w(30)=w(32) - 4 + w(30)
      w(30)=w(11)*w(30)
      w(32)=w(29) - 1
      w(33)= - u2*w(32)
      w(33)= - 2 + w(33)
      w(33)=w(2)*w(10)*w(33)
      w(30)=w(33) + w(30)
      w(30)=w(2)*w(30)
      w(33)=w(20) - 4
      w(34)= - u2*w(20)
      w(34)= - 2*w(33) + w(34)
      w(34)=u2*w(34)
      w(34)=w(34) - 12 - w(20)
      w(34)=u2*w(34)
      w(34)=4 + w(34)
      w(34)=w(3)*w(34)
      w(35)= - 8 + u2
      w(35)=u2*w(35)
      w(35)=5 + w(35)
      w(35)=u2*w(35)
      w(35)= - 2 + w(35)
      w(35)=w(5)*w(35)
      w(34)=w(35) + w(34)
      w(34)=w(12)*w(34)
      w(21)=w(34) + w(27) + w(30) + w(21)
      w(27)=w(19) - 2
      w(30)=w(20) + w(27)
      w(30)=w(30)*w(29)
      w(30)=w(30) - 4*w(27) - 5*w(20)
      w(30)=u2*w(30)
      w(19)=w(30) + w(19) + w(33)
      w(19)=u2*w(19)
      w(19)=w(19) + 8 + w(20)
      w(19)=u2*w(19)
      w(24)=w(24) - 1
      w(30)= - w(24)*w(29)
      w(28)=w(30) - 5 + w(28)
      w(28)=u2*w(28)
      w(25)=w(25) - 1
      w(28)=w(28) - w(25)
      w(28)=u2*w(28)
      w(30)=w(20) - 5
      w(28)=w(28) + w(30)
      w(28)=u2*w(28)
      w(29)= - w(31)*w(29)
      w(29)= - 3 + w(29)
      w(29)=u2*w(29)
      w(29)=4 + w(29)
      w(29)=u2*w(29)
      w(29)= - 2 + w(29)
      w(29)=w(2)*w(29)
      w(28)=w(29) + 2 + w(28)
      w(28)=w(2)*w(28)
      w(20)=w(20) + 1
      w(20)=w(20)*u2
      w(25)= - w(20) + w(25)
      w(25)=u2*w(25)
      w(25)=w(25) - w(30)
      w(25)=u2*w(25)
      w(29)=u2 - 4
      w(29)=w(29)*u2
      w(30)= - 8 - w(29)
      w(30)=u2*w(30)
      w(30)=4 + w(30)
      w(30)=w(2)*w(30)
      w(25)=w(30) - 2 + w(25)
      w(25)=w(3)*w(25)
      w(29)= - 2 - w(29)
      w(29)=w(5)*w(29)
      w(19)=w(25) + w(28) + 2*w(29) - 4 + w(19)
      w(19)=w(17)*w(19)
      w(18)=w(16)*w(32)*w(18)
      w(20)=w(20) - w(24)
      w(20)=u2*w(20)
      w(20)=1 + w(20)
      w(20)=w(16)*w(20)
      w(24)=w(1)*w(2)*w(14)
      w(20)=w(20) + w(24)
      w(20)=w(4)*w(20)
      w(23)= - w(15)*w(23)
      w(18)=w(20) + w(23) + w(19) + w(18)
      w(18)=2*w(18)
      w(19)= - im*w(26)
      w(20)=w(2)*im
      w(23)= - Pi + 2*im
      w(23)=2*w(23) + w(20)
      w(23)=w(2)*w(23)
      w(19)=w(19) + w(23)
      w(23)= - im*w(27)
      w(20)=w(20) - im - 2*Pi
      w(20)=w(2)*w(20)
      w(20)=w(20) + Pi + w(23)
      w(20)=4*w(1)*w(20)


      omega_ggg_ppp_1l_2a_mt0 =  0

      omega_ggg_pmp_1l_2a_mt0 = w(21)

      omega_qag_mpp_1l_2a_mt0 = w(19)

      omega_ggg_ppp_1l_2a_mt2 = w(22)

      omega_ggg_pmp_1l_2a_mt2 = w(18)

      omega_qag_mpp_1l_2a_mt2 = w(20)
