!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function ttbqqbrqpm(k1,k2,k3,k4,k5,k6,k7)
      implicit none
      include 'types.f'
      complex(dp):: ttbqqbrqpm
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer:: k1,k2,k3,k4,k5,k6,k7
      real(dp):: s3459,s6789,mtsq
      s3459=s(k4,k3)+s(k5,k3)
      s6789=s(k3,k6)+s(k3,k7)
      mtsq=mt**2
      ttbqqbrqpm =  + xn**(-1)*s6789**(-1) * (  - 1/(za(k1,k2))/(zb(k1,
     &    k2))/(zb(k5,k3))*za(k1,k5)*za(k6,k3)*za(k6,k7)*zb(k2,k6)*zb(
     &    k4,k5)*zb(k5,k6) - 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(
     &    k1,k5)*za(k6,k7)*za(k7,k3)*zb(k2,k7)*zb(k4,k5)*zb(k5,k6) )
      ttbqqbrqpm = ttbqqbrqpm + xn**(-1)*s3459**(-1) * ( 1/(za(k1,k2))
     &    /(zb(k1,k2))*za(k1,k3)*za(k5,k3)*za(k6,k7)*zb(k2,k6)*zb(k4,k5
     &    ) - 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(k1,k4)*za(k5,k3)
     &    *za(k6,k7)*zb(k2,k6)*zb(k4,k5)**2 )
      ttbqqbrqpm = ttbqqbrqpm + xn**(-1)*mtsq*s6789**(-1) * ( 1/(za(k1,
     &    k2))/(zb(k1,k2))*za(k1,k3)*za(k7,k3)*zb(k2,k4) + 1/(za(k1,k2)
     &    )/(zb(k1,k2))/(zb(k5,k3))*za(k1,k3)*za(k6,k7)*zb(k2,k4)*zb(k5
     &    ,k6) + 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(k1,k5)*za(k7,
     &    k3)*zb(k2,k5)*zb(k4,k5) + 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3
     &    ))*za(k1,k6)*za(k7,k3)*zb(k2,k4)*zb(k5,k6) + 1/(za(k1,k2))/(
     &    zb(k1,k2))/(zb(k5,k3))*za(k1,k7)*za(k7,k3)*zb(k2,k4)*zb(k5,k7
     &    ) )
      ttbqqbrqpm = ttbqqbrqpm + xn**(-1)*mtsq*s3459**(-1) * ( 1/(za(k1,
     &    k2))/(zb(k1,k2))/(zb(k5,k3))*za(k1,k3)*za(k6,k7)*zb(k2,k6)*
     &    zb(k4,k5) + 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(k1,k7)*
     &    za(k4,k3)*zb(k2,k4)*zb(k4,k5) )

      return
      end
