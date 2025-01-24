!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      w(1)=HZ1(1)
      w(2)=GYZ1(0)
      w(3)=GYZ1(3)
      w(4)=HZ2(0,1)
      w(5)=GYZ2(0,2)
      w(6)=GYZ2(1,0)
      w(7)=GYZ2(2,0)
      w(8)=GYZ2(3,2)
      w(9)=OneMz**(-1)
      w(10)=HZ1(0)
      w(11)=y**(-1)
      w(12)=GYZ1(2)
      w(13)=yPz**(-1)
      w(14)=w(3)*w(1)
      w(14)=w(14) - w(8)
      w(15)=w(2)*w(1)
      w(15)= - w(7) - w(4) + w(15) - w(14) + w(6) - w(5)
      w(16)=CF*w(15)
      w(17)=w(10) + w(1)
      w(18)=w(10) + 3._dp/2._dp
      w(18)=w(18)*w(12)
      w(17)= - w(18) - w(14) + 3._dp/2._dp*w(17)
      w(18)=1 + 2*w(10)
      w(18)=w(18)*w(12)
      w(18)=w(1) - w(18) - 2*w(14)
      w(19)=w(18) + w(10)
      w(19)=w(19)*z
      w(20)=w(12)*w(10)
      w(14)=w(14) + w(20)
      w(20)= - w(11)*w(14)*z**2
      w(21)=w(9)*w(10)
      w(22)=w(21) - w(10)
      w(20)=w(20) + w(19) - w(22)
      w(20)=w(11)*w(20)
      w(21)=1._dp/2._dp*w(21)
      w(23)= - w(21) - 1._dp/2._dp - w(10)
      w(23)=w(9)*w(23)
      w(20)=w(20) + w(23) - 7._dp/2._dp + w(17)
      w(20)=CF*w(20)
      w(23)=w(12) - w(1)
      w(24)=w(23)*w(13)
      w(25)=w(24) + 1
      w(26)=z*w(13)
      w(26)=1._dp/2._dp*w(26)
      w(25)=w(25)*w(26)
      w(15)= - w(25) - w(24) - 4 + w(15)
      w(15)=CF*w(15)
      w(14)=w(14)*w(11)*z
      w(26)= - z + 1
      w(26)=w(26)*w(14)
      w(18)=w(26) + w(19) - w(18)
      w(18)=w(11)*w(18)
      w(17)=w(18) - w(25) + w(21) + 2*w(24) + 1 + w(17)
      w(17)=CF*w(17)
      w(18)=CF*w(25)
      w(14)=w(14) + w(23) + w(22)
      w(14)=w(11)*w(14)
      w(19)=1._dp/2._dp*w(9)
      w(21)=1 + w(22)
      w(19)=w(21)*w(19)
      w(14)=w(14) + w(25) + w(19) - 1._dp/2._dp - w(24)
      w(14)=CF*w(14)


      alphaG1WL2= 0

      alphaG2WL2= 0

      alphaG4WL2= 0

      betaG1WL2= 0

      betaG2WL2= 0

      betaG4WL2= 0

      gammaG1WL2= 0

      gammaG2WL2= 0

      gammaG4WL2= 0

      alphaG1WL1=w(16)

      alphaG2WL1=w(20)

      betaG1WL1=w(15)

      betaG2WL1=w(17)

      gammaG1WL1=w(18)

      gammaG2WL1=w(14)

      alphaG1WL0= 0

      alphaG2WL0=1

      betaG1WL0=1

      betaG2WL0= 0

      gammaG1WL0= 0

      gammaG2WL0= 0
