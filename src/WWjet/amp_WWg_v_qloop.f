!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine amp_WWg_v_qloop(j1,j2,j3,j4,j5,j6,j7,
     & ampq,ampq12,ampqflip,ampq12flip,
     & ampax,ampax12,ampaxflip,ampax12flip)
c--- master routine for assembling amplitude for
c---  ub(i1) + u(i2) -> nu(i3) + e+(i4) + e-(i5) + nu~(i6) + g(i7)
c---
c--- amplitude is returned as amp1a(h1,h2) and amp1b(h1,h2) with h1,h2
c--- the helicities of the quark and gluon respectively
c--- ampa refers to diagrams with both W's attached to the quark line
c--- ampb refers to diagrams with Z->WW
c--- routine also returns the corresponding LO amplidues amp0a,amp0b
c---
c--- Flavor of quark line is passed in as 'Qid'
      implicit none
      include 'types.f'
      integer j1,j2,j3,j4,j5,j6,j7,polq
      complex(dp):: ampq(2,2),ampq12(2,2),ampqflip(2,2),ampq12flip(2,2),
     & ampax(2,2),ampax12(2,2),ampaxflip(2,2),ampax12flip(2,2)

      call AWWjetq(j1,j2,j3,j4,j5,j6,j7,ampq,ampax)

c Amplitudes for all relevant crossings can now be related to these
      do polq=1,2
      ampq12(polq,:)=ampq(3-polq,:)
      ampax12(polq,:)=ampax(3-polq,:)
      ampqflip(polq,:)=ampq(polq,:)
      ampaxflip(polq,:)=-ampax(polq,:)
      ampq12flip(polq,:)=ampq(3-polq,:)
      ampax12flip(polq,:)=-ampax(3-polq,:)
      enddo

      return
      end

