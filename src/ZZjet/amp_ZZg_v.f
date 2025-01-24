!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine amp_ZZg_v(p,j1,j2,j3,j4,j5,j6,j7,amp0a,amp0sr,amp1a,amp1ax56,amp1ax34,amp1sr)
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
      include 'srdiags.f'
      integer:: j1,j2,j3,j4,j5,j6,j7,h34,h56
      real(dp):: p(mxpart,4)
      complex(dp):: amp0a(2,2,2,2),amp1a(2,2,2,2),
     & amp1ax56(2,2,2,2),amp1ax34(2,2,2,2),amp0sr(2,2,2,2,2),amp1sr(2,2,2,2,2),tmp(2,2,2,2)
c--- a and q amplitudes no longer computed together, for efficiency
c      call AZZjetaq(p,j1,j2,j3,j4,j5,j6,j7,amp0a,amp0sr,amp1a,amp1q,amp1sr)
c      call AZZjetq(p,j1,j2,j3,j4,j5,j6,j7,amp1q)
c--- a amplitudes only
      call AZZjeta(p,j1,j2,j3,j4,j5,j6,j7,amp0a,amp0sr,amp1a,amp1sr)

      if (srdiags) then
        call AZZjet_qlooptri_sr(j1,j2,j3,j4,j6,j5,j7,tmp) ! Note j6,j5
        amp1sr(1,:,:,:,:)=amp1sr(1,:,:,:,:)+tmp(:,:,:,:)
        amp1sr(2,:,:,:,:)=amp1sr(2,:,:,:,:)+tmp(:,:,:,:)
      endif

c--- axial triangle contributions
      call AZZjetax(j1,j2,j3,j4,j6,j5,j7,amp1ax56)     ! Note j6,j5
      call AZZjetax(j1,j2,j6,j5,j3,j4,j7,tmp)          ! Note j6,j5
      do h34=1,2
      do h56=1,2
      amp1ax34(:,h34,h56,:)=tmp(:,h56,h34,:)
      enddo
      enddo

      return
      end

