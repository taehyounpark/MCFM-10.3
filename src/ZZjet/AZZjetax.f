!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AZZjetax(j1,j2,j3,j4,j5,j6,j7,A7ax)
c--- Compute amplitude A^q_7 virtual, for both positive
c--- and negative gluon helicities, and for both lepton helicities for Z decay

c--- Result for diagrams of the form:

c        u1
c           |      l3
c           | Z   /              l6
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
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: A7ax(2,2,2,2),bit(2,2)

c labelling of A7q is (pol12, pol34, pol56, polg)

c labelling of bit is (polg, pol56); pol12, pol34 are naturally left-handed
      call AVZjet_tri(j1,j2,j3,j4,j5,j6,j7,bit)
      A7ax(1,1,:,1)=bit(1,:)
      A7ax(1,1,:,2)=bit(2,:)
      call AVZjet_tri(j1,j2,j4,j3,j5,j6,j7,bit)
      A7ax(1,2,:,1)=bit(1,:)
      A7ax(1,2,:,2)=bit(2,:)
      call AVZjet_tri(j2,j1,j3,j4,j5,j6,j7,bit)
      A7ax(2,1,:,1)=-bit(1,:)
      A7ax(2,1,:,2)=-bit(2,:)
      call AVZjet_tri(j2,j1,j4,j3,j5,j6,j7,bit)
      A7ax(2,2,:,1)=-bit(1,:)
      A7ax(2,2,:,2)=-bit(2,:)

      return
      end


