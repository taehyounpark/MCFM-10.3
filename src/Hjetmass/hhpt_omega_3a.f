!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      !complex(dp) w(25)

      w(1)=s**(-1)
      w(2)=Log(1 - u3)
      w(3)=Log(mt2*t**(-1))
      w(4)=Log(u3)
      w(5)=1/( - 1 + u3)
      w(6)=1/(s - 2*s*u3 + s*u3**2)
      w(7)=w(2)**2
      w(8)=Pi**2
      w(9)=2*w(8) + w(7)
      w(9)=w(1)*w(9)
      w(10)=im*Pi
      w(11)=w(10)*w(1)
      w(12)=2*w(2)
      w(13)=w(12)*w(11)
      w(14)=w(2)*w(1)
      w(11)=w(14) - w(11)
      w(14)=w(3)*w(1)
      w(11)=2*w(11) - 3*w(14)
      w(11)=w(3)*w(11)
      w(15)=Pi*w(1)
      w(16)=im*w(1)
      w(17)=w(2)*w(16)
      w(15)=w(15) + w(17)
      w(15)=im*w(15)
      w(14)=w(15) + w(14)
      w(15)=w(4)*w(1)
      w(14)=2*w(14) + w(15)
      w(14)=w(4)*w(14)
      w(9)=w(14) + w(11) + w(13) + w(9)
      w(9)=2*mt2*w(9)
      w(11)=w(8) + w(7) + 4
      w(13)=im + Pi
      w(13)=w(13)*im
      w(13)=w(13) + w(3)
      w(14)=2*w(4)
      w(15)= - w(14) - 3 + 4*w(2) + w(13)
      w(15)=w(4)*w(15)
      w(17)=w(2)*Pi
      w(18)=w(17) + 2*Pi
      w(19)=im*w(18)
      w(20)=2*w(3)
      w(21)= - w(20) + 10 + 3*w(2)
      w(21)=w(3)*w(21)
      w(15)=w(15) + w(21) - 2*w(11) + 3*w(19)
      w(19)=mt2*w(6)
      w(15)=w(15)*w(19)
      w(21)=im*w(2)
      w(17)= - w(21) - 3*Pi - w(17)
      w(17)=im*w(17)
      w(12)=w(12) - w(4)
      w(13)= - w(12) - w(13)
      w(13)=w(4)*w(13)
      w(21)=w(3) - 5 - w(2)
      w(21)=w(3)*w(21)
      w(11)=w(13) + w(21) + w(17) + w(11)
      w(11)=w(11)*w(19)
      w(7)=w(7)*w(5)
      w(13)=w(8)*w(5)
      w(11)=w(11) + w(13) + w(7)
      w(11)=u3*w(11)
      w(17)=w(2)*w(5)
      w(21)=w(5) - w(17)
      w(21)=w(2)*w(21)
      w(21)= - w(13) + w(21)
      w(12)= - w(4)*w(12)
      w(11)=w(11) + w(15) + 2*w(21) + w(12)
      w(11)=u3*w(11)
      w(11)=w(11) + 6*w(5) + w(13)
      w(12)=2*im
      w(13)= - w(18)*w(12)
      w(15)= - w(2) + w(10) + 2
      w(18)=4*w(3)
      w(21)= - w(4) + w(18) + w(15)
      w(21)=w(4)*w(21)
      w(22)= - 4 - w(2)
      w(22)=w(2)*w(22)
      w(15)= - w(15) - 3*w(3)
      w(15)=w(3)*w(15)
      w(8)=w(21) + w(15) + w(13) + w(22) - 4 + w(8)
      w(13)=2*w(19)
      w(8)=w(8)*w(13)
      w(15)=w(10)*w(5)
      w(21)=4*w(5)
      w(22)=w(3)*w(5)
      w(15)=3.D0/2.D0*w(22) + w(15) + w(21)
      w(23)= - w(17) + w(15)
      w(23)=w(3)*w(23)
      w(24)= - 1 + w(2)
      w(14)= - w(14) + 4*w(24) - w(22)
      w(14)=w(4)*w(14)
      w(24)=w(17)*w(10)
      w(25)= - 8*w(5) + 5.D0/2.D0*w(17)
      w(25)=w(2)*w(25)
      w(8)=w(8) + w(14) + w(23) - w(24) + w(25) + 2*w(11)
      w(8)=u3*w(8)
      w(11)= - w(17) - w(15)
      w(11)=w(3)*w(11)
      w(10)=w(10) + w(2)
      w(14)=w(10) - 3
      w(15)=w(18) + w(14)
      w(15)=w(3)*w(15)
      w(14)=4*w(4) - 8*w(3) - w(14)
      w(14)=w(4)*w(14)
      w(14)=w(14) + w(15) + 8 + w(10)
      w(13)=w(14)*w(13)
      w(10)=1.D0/2.D0*w(4) + 3*w(22) - w(10)
      w(10)=w(4)*w(10)
      w(14)=12*w(5)
      w(7)=w(8) + w(13) + w(10) + w(11) + w(24) - w(14) + 1.D0/2.D0*
     & w(7)
      w(7)=u3*w(7)
      w(8)=w(21) + w(22)
      w(8)=w(3)*w(8)
      w(10)=w(3) - 1
      w(10)=w(10)*w(3)
      w(11)= - w(4) + w(20) - 1
      w(11)=w(11)*w(4)
      w(10)= - w(11) + w(10) + 2
      w(11)= - w(10)*w(19)
      w(13)=2 - w(22)
      w(13)=2*w(13) - w(4)
      w(13)=w(4)*w(13)
      w(7)=w(7) + 4*w(11) + w(13) + w(14) + w(8)
      w(8)=w(3)*im
      w(11)=4 + w(3)
      w(11)=w(11)*w(8)
      w(8)= - w(12) - w(8)
      w(12)=w(4)*im
      w(8)=2*w(8) + w(12)
      w(8)=w(4)*w(8)
      w(10)=mt2*w(16)*w(10)
      w(8)=4*w(10) + w(8) + 12*im + w(11)


      omega_ggg_ppp_1l_3a = w(9)

      omega_ggg_pmp_1l_3a = w(7)

      omega_qag_mpp_1l_3a = w(8)
