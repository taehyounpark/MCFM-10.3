!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function D3six(k1,k2,k3,k4,k5,k6)
      use loopI3_generic
      use loopI4_generic
      implicit none
      include 'types.f'
      complex(dp):: D3six

c----- six dimensional box corresponding to
c----- loopI4(s56,0._dp,s34,0._dp,s156,s134,0._dp,mtsq,mtsq,0._dp,musq,e)
c----- multiplied by a factor of -C(0)/2
      include 'mxpart.f'
      include 'masses.f'
      include 'scale.f'
      include 'sprods_com.f'
      include 'scalarselect.f'
      integer:: k1,k2,k3,k4,k5,k6
      complex(dp):: IntC(4),IntD
      real(dp):: s12,s34,s56,s134,s156,C(0:4),mtsq,Delta
      integer:: e
      mtsq=mt**2
      s12=s(k1,k2)
      s34=s(k3,k4)
      s56=s(k5,k6)
      s134=s(k1,k3)+s(k1,k4)+s(k3,k4)
      s156=s(k1,k5)+s(k1,k6)+s(k5,k6)
      Delta=s12*mtsq-s34*s56+s134*s156
      C(1)=2._dp*(s134-s34)*(2._dp*s12*mtsq-Delta)/Delta**2
      C(2)=2._dp*(s34-s156)/Delta
      C(3)=2._dp*(s56-s134)/Delta
      C(4)=2._dp*(s156-s56)*(2._dp*s12*mtsq-Delta)/Delta**2
      C(0)=4._dp*s12*(s134*s156-s34*s56)/Delta**2
      e=0
      IntC(1)=loopI3(s134,0._dp,s34,0._dp,mtsq,mtsq,musq,e)    !C7
      IntC(2)=loopI3(0._dp,s34,s156,0._dp,0._dp,mtsq,musq,e)     !C5
      IntC(3)=loopI3(0._dp,s56,s134,0._dp,0._dp,mtsq,musq,e)     !C3
      IntC(4)=loopI3(s156,0._dp,s56,0._dp,mtsq,mtsq,musq,e)    !C9
      IntD=loopI4(s56,0._dp,s34,0._dp,s156,s134,0._dp,mtsq,mtsq,0._dp,musq,e)
      D3six=0.5_dp*(C(1)*IntC(1)+C(2)*IntC(2)+C(3)*IntC(3)+C(4)*IntC(4)
     & +2._dp*IntD)

c---debug
c---real six dim box
c      D3six=-(C(1)*IntC(1)+C(2)*IntC(2)+C(3)*IntC(3)+C(4)*IntC(4)
c     & +2._dp*IntD)/C(0)
c---debug
      return
      end

