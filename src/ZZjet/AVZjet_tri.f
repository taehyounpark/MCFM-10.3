!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AVZjet_tri(j1,j2,j3,j4,j5,j6,j7,A7q)
c--- Compute amplitude A^q_7 virtual, for both positive
c--- and negative gluon helicities, and for both lepton helicities for Z decay

c--- Result for diagrams of the form:

c        u1
c           |      l3
c           |  V  /              l6
c           |~~~~~           Z  /
c           |     \         ~~~~
c           |      a4     /|    \
c           ^            / |     a5
c           |           /  |
c           |ooooooooooo   |
c           |           \  |
c           |            \ |
c           |             \|
c        d2                 ooooo g7

c--- Obtained by contracting result in Appendix F of 1610.02189
c--- with the appropriate currents

      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'scale.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer j1,j2,j3,j4,j5,j6,j7
      real(dp):: s1234
      complex(dp):: A7q(2,2),F1anom,G1,amp_qloop_ax56_mmmp

      s1234=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)

      G1=F1anom(s1234,s(j5,j6),mt**2,musq)-F1anom(s1234,s(j5,j6),mb**2,musq)

c Indices label hg, h56 (with h12=h34=1)
      A7q(2,1)=+amp_qloop_ax56_mmmp(j1,j2,j3,j4,j5,j6,j7,za,zb)
      A7q(2,2)=+amp_qloop_ax56_mmmp(j1,j2,j3,j4,j6,j5,j7,za,zb)
      A7q(1,1)=-amp_qloop_ax56_mmmp(j2,j1,j4,j3,j6,j5,j7,zb,za)
      A7q(1,2)=-amp_qloop_ax56_mmmp(j2,j1,j4,j3,j5,j6,j7,zb,za)
      A7q(:,:)=A7q(:,:)*G1

c--- overall factor remaining: gsq*e/16/pisq delta(A,B)
c---                          *gwsq*e/2*vLR tB
c---                          = ason2pi * 2*gwsq*esq  * vLR * tA

      return
      end


      function amp_qloop_ax56_mmmp(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zcouple_cms.f'
      complex(dp):: amp_qloop_ax56_mmmp,zab2
      real(dp):: s134,s234
      integer p1,p2,p3,p4,p5,p6,p7

c--- statement functions
      zab2(p1,p3,p4,p2)=za(p1,p3)*zb(p3,p2)+za(p1,p4)*zb(p4,p2)

      s134=s(p1,p3)+s(p1,p4)+s(p3,p4)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)

      amp_qloop_ax56_mmmp = zb(p7,p6)/(s(p3,p4)*s(p5,p6))
     & *(zab2(p5,p1,p3,p4)*zb(p2,p7)*za(p1,p3)/s134
     &  -za(p5,p1)*zab2(p3,p2,p4,p7)*zb(p4,p2)/s234)

      amp_qloop_ax56_mmmp = amp_qloop_ax56_mmmp/(sqrt(zxw*(cone-zxw)))


      return
      end


