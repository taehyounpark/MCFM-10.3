!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine GLYF1(qsq,F1)
c     Elaborated using the form results of 1403.6451v2
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
      real(dp):: F1(0:1,0:2)
c      F1,F2 index one transition flavor gg=0,qq=1
c      F1,F2 index2,power of Lb
      real(dp)::qsq,Lqsq,F11(0:1)
      real(dp),parameter::C(0:1)=[CA,CF]
      integer:: i

c     initialize to zero
      F1(:,:)=0

c     A bit of a question here whether Lqsq is real or complex
      Lqsq=log(qsq/musq)
      do i=0,1
      F11(i)=4*C(i)

      F1(i,2)=-F11(i)
      F1(i,1)=-F11(i)*Lqsq
      F1(i,0)=0

      enddo
      return
      end

