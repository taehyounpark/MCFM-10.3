!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine amp_WZg_v(p,j1,j2,j3,j4,j5,j6,j7,amp0ad,amp0au,amp0b,amp0sr,amp1ad,amp1au,amp1b,amp1q,amp1sr)
c--- master routine for assembling amplitude for
c---  u(i1) + dbar(i2) -> e-(i3) + nu~(i4) + mu+(i5) + mu-(i6) + g(i7)
c---
c--- amplitude is returned as amp1a(h1) with h1 the helicity of the gluon
c--- ampa refers to diagrams with W and Z attached to the quark line,
c---  u,d requires dressing with vL(u,d) respectively (quark is always LH)
c--- ampb refers to diagrams with W->WZ
c--- routine also returns the corresponding LO amplidues amp0a,amp0b
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'srdiags.f'
      integer j1,j2,j3,j4,j5,j6,j7
      real(dp):: p(mxpart,4)
      complex(dp)::
     & amp0b(2,2),amp1b(2,2),
     & amp0ad(2,2),amp0au(2,2),amp1ad(2,2),amp1au(2,2),amp1q(2,2),
     & amp0sr(2,2),amp1sr(2,2)

      amp0ad(:,:)=czip
      amp0au(:,:)=czip
      amp0b(:,:)=czip
      amp1ad(:,:)=czip
      amp1au(:,:)=czip
      amp1b(:,:)=czip

c--- a amplitudes
      call AWmZjeta(p,j1,j2,j3,j4,j5,j6,j7,amp0ad,amp0au,amp1ad,amp1au)

c--- b amplitude: both left-handed and right-handed quarks
      call AWmZjetb(j1,j2,j3,j4,j5,j6,j7,amp0b,amp1b)

      if (srdiags) then
        call AWmZjetsr(j1,j2,j3,j4,j5,j6,j7,amp0sr,amp1sr)
      else
        amp0sr(:,:)=czip
        amp1sr(:,:)=czip
      endif

c--- axial triangle contribution
      call AVZjet_tri(j1,j2,j3,j4,j6,j5,j7,amp1q)     ! Note j6,j5

      return
      end

