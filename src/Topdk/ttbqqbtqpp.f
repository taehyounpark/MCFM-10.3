!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function ttbqqbtqpp(k1,k2,k3,k4,k5,k6,k7)
      implicit none
      include 'types.f'
      complex(dp):: ttbqqbtqpp
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer:: k1,k2,k3,k4,k5,k6,k7
      real(dp):: s129,s6789,mtsq
      s129=s(k1,k2)+s(k1,k3)+s(k2,k3)
      s6789=s(k3,k6)+s(k3,k7)
      mtsq=mt**2
      ttbqqbtqpp =  + s129**(-1) * (  - 1/(za(k1,k2))/(za(k5,k3))/(zb(
     &    k1,k2))*za(k1,k5)**2*za(k6,k7)*zb(k1,k3)*zb(k2,k6)*zb(k4,k5)
     &     + 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k5)**2*za(k6,
     &    k7)*zb(k1,k6)*zb(k2,k3)*zb(k4,k5) - 1/(za(k2,k3))*za(k1,k5)*
     &    za(k6,k7)*zb(k4,k5)*zb(k6,k3) - 1/(za(k2,k3))/(za(k5,k3))*za(
     &    k1,k5)*za(k2,k5)*za(k6,k7)*zb(k2,k6)*zb(k4,k5) )
      ttbqqbtqpp = ttbqqbtqpp + s6789**(-1) * (  - 1/(za(k1,k2))/(za(k5
     &    ,k3))/(zb(k1,k2))*za(k1,k5)*za(k5,k6)*za(k6,k7)*zb(k2,k6)*zb(
     &    k4,k5)*zb(k6,k3) - 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(
     &    k1,k5)*za(k5,k7)*za(k6,k7)*zb(k2,k7)*zb(k4,k5)*zb(k6,k3) -
     &    1/(za(k1,k2))/(zb(k1,k2))*za(k1,k5)*za(k6,k7)*zb(k2,k3)*zb(k4
     &    ,k5)*zb(k6,k3) )
      ttbqqbtqpp = ttbqqbtqpp + mtsq*s129**(-1) * ( 1/(za(k1,k2))/(za(
     &    k5,k3))/(zb(k1,k2))*za(k1,k3)*za(k5,k7)*zb(k2,k3)*zb(k4,k3)
     &     + 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k5)*za(k1,k7)*
     &    zb(k1,k3)*zb(k2,k4) - 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*
     &    za(k1,k5)*za(k1,k7)*zb(k1,k4)*zb(k2,k3) - 1/(za(k1,k2))/(za(
     &    k5,k3))/(zb(k1,k2))*za(k1,k5)*za(k2,k7)*zb(k2,k3)*zb(k2,k4)
     &     + 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k7)*za(k2,k5)*
     &    zb(k2,k3)*zb(k2,k4) + 1/(za(k2,k3))*za(k1,k7)*zb(k4,k3) + 1/(
     &    za(k2,k3))/(za(k5,k3))*za(k1,k7)*za(k2,k5)*zb(k2,k4) )
      ttbqqbtqpp = ttbqqbtqpp + mtsq*s6789**(-1) * (  - 1/(za(k1,k2))/(
     &    za(k5,k3))/(zb(k1,k2))*za(k1,k5)*za(k5,k7)*zb(k2,k3)*zb(k4,k5
     &    ) - 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k5)*za(k6,k7)
     &    *zb(k2,k4)*zb(k6,k3) + 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*
     &    za(k1,k6)*za(k5,k7)*zb(k2,k4)*zb(k6,k3) + 1/(za(k1,k2))/(za(
     &    k5,k3))/(zb(k1,k2))*za(k1,k7)*za(k5,k7)*zb(k2,k4)*zb(k7,k3) )

      return
      end
