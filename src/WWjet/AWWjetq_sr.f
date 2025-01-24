!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWWjetq_sr(j1,j2,j3,j4,j5,j6,j7,A7ax)
c--- closed fermion loop contributions for singly-resonant diagrams
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: A7ax(2,2),AWWjet_qlooptri_sr

c--- positive helicity gluon
      A7ax(1,2)=-AWWjet_qlooptri_sr(j1,j2,j6,j5,j3,j4,j7,za,zb)

c--- flip helicity of quark line by interchanging j1 and j2
      A7ax(2,2)=-AWWjet_qlooptri_sr(j2,j1,j6,j5,j3,j4,j7,za,zb)

c--- negative helicity gluon
      A7ax(1,1)=AWWjet_qlooptri_sr(j2,j1,j4,j3,j5,j6,j7,zb,za)

c--- flip helicity of quark line by interchanging j1 and j2
      A7ax(2,1)=AWWjet_qlooptri_sr(j1,j2,j4,j3,j5,j6,j7,zb,za)

      return
      end


