!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWWjetq(j1,j2,j3,j4,j5,j6,j7,A7q,A7ax)
c--- closed fermion loop contributions
c--- includes nf (=5) flavors of massless quark box diagrams
c--- and one (t,b) generation of massive triangle diagrams
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'WWjetlabels.f'
      include 'zprods_com.f'
      include 'verbose.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: A7q(2,2),A7ax(2,2),AWWjet_qlooptri

      call WWjetcomputescalars(j1,j2,j3,j4,j5,j6,j7,scints)

      verbose=.false.

c--- positive helicity gluon
      call AWWjet_qloop_new('mp',.false.,j1,j2,j3,j4,j6,j5,j7,za,zb,scints,A7q(1,2))
      A7ax(1,2)=-AWWjet_qlooptri(j1,j2,j6,j5,j3,j4,j7,za,zb)

c--- flip helicity of quark line by interchanging j1 and j2
      call AWWjet_qloop_new('pp',.false.,j2,j1,j3,j4,j6,j5,j7,za,zb,scints,A7q(2,2))
      A7ax(2,2)=-AWWjet_qlooptri(j2,j1,j6,j5,j3,j4,j7,za,zb)

c--- negative helicity gluon
      call AWWjet_qloop_new('mm',.true.,j2,j1,j5,j6,j4,j3,j7,zb,za,scints,A7q(1,1))
      A7ax(1,1)=AWWjet_qlooptri(j2,j1,j4,j3,j5,j6,j7,zb,za)

c--- flip helicity of quark line by interchanging j1 and j2
      call AWWjet_qloop_new('pm',.true.,j1,j2,j5,j6,j4,j3,j7,zb,za,scints,A7q(2,1))
      A7ax(2,1)=AWWjet_qlooptri(j1,j2,j4,j3,j5,j6,j7,zb,za)

c--- hard-coded for two massless generations of quarks
      A7q(:,:)=A7q(:,:)*4._dp

      return
      end


