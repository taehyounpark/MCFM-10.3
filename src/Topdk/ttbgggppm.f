!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function ttbgggppm(i1,i2,i3,i4,i5,i6,i7)
      implicit none
      include 'types.f'
      complex(dp):: ttbgggppm
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      real(dp):: s129,s1345,s6789,mtsq
      s129=s(i1,i2)+s(i1,i3)+s(i2,i3)
      s1345=s(i1,i4)+s(i1,i5)
      s6789=s(i3,i6)+s(i3,i7)
      mtsq=mt**2
      ttbgggppm =  + s129**(-1) * (  - 1/(za(i1,i2))/(za(i1,i3))/(zb(i6
     &    ,i3))*za(i2,i3)*za(i5,i3)*za(i6,i7)*zb(i2,i6)**2*zb(i4,i5) -
     &    1/(za(i1,i2))/(za(i2,i3))/(zb(i6,i3))*za(i1,i3)*za(i5,i3)*za(
     &    i6,i7)*zb(i1,i6)**2*zb(i4,i5) - 2._dp/(za(i1,i2))/(zb(i6,i3))*
     &    za(i5,i3)*za(i6,i7)*zb(i1,i6)*zb(i2,i6)*zb(i4,i5) - 1/(za(i1,
     &    i3))/(zb(i2,i3))/(zb(i6,i3))*za(i5,i3)*za(i6,i7)*zb(i1,i2)*
     &    zb(i2,i6)**2*zb(i4,i5) - 1/(za(i2,i3))/(zb(i2,i3))/(zb(i6,i3)
     &    )*za(i5,i3)*za(i6,i7)*zb(i1,i2)*zb(i1,i6)*zb(i2,i6)*zb(i4,i5)
     &     )
      ttbgggppm = ttbgggppm + s1345**(-1) * ( 1/(za(i1,i3))/(za(i2,i3))
     &    /(zb(i2,i3))/(zb(i6,i3))*za(i4,i3)*za(i5,i3)*za(i6,i7)*zb(i1,
     &    i4)*zb(i2,i6)**2*zb(i4,i5) + 1/(za(i1,i3))/(za(i2,i3))/(zb(i2
     &    ,i3))/(zb(i6,i3))*za(i5,i3)**2*za(i6,i7)*zb(i1,i5)*zb(i2,i6)
     &    **2*zb(i4,i5) )
      ttbgggppm = ttbgggppm + mtsq*s129**(-1) * ( 1/(za(i1,i2))/(za(i1,
     &    i3))/(zb(i6,i3))*za(i2,i3)*za(i7,i3)*zb(i2,i4)*zb(i2,i6) +
     &    1/(za(i1,i2))/(za(i2,i3))/(zb(i6,i3))*za(i1,i3)*za(i7,i3)*zb(
     &    i1,i4)*zb(i1,i6) + 1/(za(i1,i2))/(zb(i6,i3))*za(i7,i3)*zb(i1,
     &    i4)*zb(i2,i6) + 1/(za(i1,i2))/(zb(i6,i3))*za(i7,i3)*zb(i1,i6)
     &    *zb(i2,i4) + 1/(za(i1,i3))/(zb(i2,i3))/(zb(i6,i3))*za(i7,i3)*
     &    zb(i1,i2)*zb(i2,i4)*zb(i2,i6) + 1/(za(i2,i3))/(zb(i2,i3))/(
     &    zb(i6,i3))*za(i7,i3)*zb(i1,i2)*zb(i1,i4)*zb(i2,i6) )
      ttbgggppm = ttbgggppm + mtsq*s1345**(-1)*s6789**(-1) * ( 1/(za(i1
     &    ,i3))/(za(i2,i3))/(zb(i6,i3))*za(i4,i3)*za(i5,i3)*za(i7,i3)*
     &    zb(i1,i4)*zb(i2,i6)*zb(i4,i5) + 1/(za(i1,i3))/(za(i2,i3))/(
     &    zb(i6,i3))*za(i4,i3)*za(i7,i3)**2*zb(i1,i4)*zb(i2,i4)*zb(i6,
     &    i7) + 1/(za(i1,i3))/(za(i2,i3))/(zb(i6,i3))*za(i5,i3)**2*za(
     &    i7,i3)*zb(i1,i5)*zb(i2,i6)*zb(i4,i5) - 1/(za(i1,i3))/(za(i2,
     &    i3))/(zb(i6,i3))*za(i5,i3)*za(i7,i3)**2*zb(i1,i2)*zb(i4,i5)*
     &    zb(i6,i7) + 1/(za(i1,i3))/(za(i2,i3))/(zb(i6,i3))*za(i5,i3)*
     &    za(i7,i3)**2*zb(i1,i4)*zb(i2,i5)*zb(i6,i7) - 1/(za(i2,i3))/(
     &    zb(i6,i3))*za(i7,i3)**2*zb(i1,i2)*zb(i1,i4)*zb(i6,i7) )
      ttbgggppm = ttbgggppm + mtsq*s1345**(-1) * (  - 1/(za(i1,i3))/(
     &    za(i2,i3))/(zb(i2,i3))/(zb(i6,i3))*za(i4,i3)*za(i7,i3)*zb(i1,
     &    i4)*zb(i2,i4)*zb(i2,i6) + 1/(za(i1,i3))/(za(i2,i3))/(zb(i2,i3
     &    ))/(zb(i6,i3))*za(i5,i3)*za(i7,i3)*zb(i1,i2)*zb(i2,i6)*zb(i4,
     &    i5) - 1/(za(i1,i3))/(za(i2,i3))/(zb(i2,i3))/(zb(i6,i3))*za(i5
     &    ,i3)*za(i7,i3)*zb(i1,i4)*zb(i2,i5)*zb(i2,i6) + 1/(za(i2,i3))
     &    /(zb(i2,i3))/(zb(i6,i3))*za(i7,i3)*zb(i1,i2)*zb(i1,i4)*zb(i2,
     &    i6) )
      ttbgggppm = ttbgggppm + mtsq*s6789**(-1) * ( 1/(za(i1,i2))/(za(i1
     &    ,i3))/(zb(i6,i3))*za(i5,i3)*za(i7,i3)*zb(i2,i6)*zb(i4,i5) +
     &    1/(za(i1,i2))/(za(i1,i3))/(zb(i6,i3))*za(i7,i3)**2*zb(i2,i4)*
     &    zb(i6,i7) + 1/(za(i1,i2))/(za(i2,i3))/(zb(i6,i3))*za(i5,i3)*
     &    za(i7,i3)*zb(i1,i6)*zb(i4,i5) + 1/(za(i1,i2))/(za(i2,i3))/(
     &    zb(i6,i3))*za(i7,i3)**2*zb(i1,i4)*zb(i6,i7) )

      return
      end
