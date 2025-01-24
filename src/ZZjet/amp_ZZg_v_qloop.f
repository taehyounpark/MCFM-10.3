!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine amp_ZZg_v_qloop(p,j1,j2,j3,j4,j5,j6,j7,
     & amp1q,amp1qmLL,amp1qmLR,amp1q12,amp1qmLL12,amp1qmLR12)
c--- master routine for assembling amplitude for
c---  u(i1) + ubar(i2) -> e-(i3) + e+(i4) + mu+(i5) + mu-(i6) + g(i7)
c---
c--- amplitude is returned as amp1a(h1) with h1 the helicity of the gluon
c--- ampa refers to diagrams with both Z's attached to the quark line,
c--- ampq refers to diagrams with Z's attached to (box) loop of massless quarks
c--- ampax56 is diagrams with Z(56) in axial anom. triangle, Z(34) attached to quarks
c--- ampax34 is diagrams with Z(34) in axial anom. triangle, Z(56) attached to quarks
c--- routine also returns the corresponding LO amplidues amp0a
      implicit none
      include 'types.f'
      include 'mxpart.f'
      integer j1,j2,j3,j4,j5,j6,j7,hq
      real(dp):: p(mxpart,4)
      complex(dp):: amp1q(2,2,2,2),amp1q12(2,2,2,2)
      complex(dp):: amp1qmLL(2,2,2,2),amp1qmLL12(2,2,2,2)
      complex(dp):: amp1qmLR(2,2,2,2),amp1qmLR12(2,2,2,2)

c--- q amplitudes
      call AZZjetq(p,j1,j2,j3,j4,j5,j6,j7,amp1q,amp1qmLL,amp1qmLR)

c amplitudes with 1 and 2 swapped
      do hq=1,2
      amp1q12(hq,:,:,:)=amp1q(3-hq,:,:,:)
      amp1qmLL12(hq,:,:,:)=amp1qmLL(3-hq,:,:,:)
      amp1qmLR12(hq,:,:,:)=amp1qmLR(3-hq,:,:,:)
      enddo

      return
      end

