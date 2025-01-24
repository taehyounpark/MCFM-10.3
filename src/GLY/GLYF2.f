!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine GLYF2(qsq,F2)
c     Elaborated using the form results of 1403.6451v2
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nfl.f'
      include 'zeta.f'
      include 'scale.f'
      real(dp):: F2(0:1,0:4)
c      F2 index one transition flavor gg=0,qq=1
c      F2 index2,power of Lb

      real(dp)::qsq,Lqsq,F11(0:1),F20(0:1),F21(0:1),F22(0:1)
      real(dp),parameter::C(0:1)=[CA,CF]
      integer:: i

c     initialize to zero
      F2(:,:)=0

c     A bit of a question here whether Lqsq is real or complex
      Lqsq=log(qsq/musq)
      do i=0,1
      F11(i)=4*C(i)
      F22(i)=C(i)*(22/3._dp*CA-8/3._dp*nfl*TF)
      F21(i)=C(i)*(268/9._dp*CA-80/9._dp*nfl*TF-8._dp*CA*zeta2)
      F20(i)=C(i)*(808/27._dp*CA-224/27._dp*nfl*TF-28*CA*zeta3)

      F2(i,4)=0.5_dp*F11(i)**2
      F2(i,3)=F11(i)**2*Lqsq-F22(i)
      F2(i,2)=0.5_dp*F11(i)**2*Lqsq**2-F21(i)-F22(i)*Lqsq
      F2(i,1)=-F20(i)-F21(i)*Lqsq
      F2(i,0)=-F20(i)*Lqsq
      enddo
      return
      end

