!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine couplz_cms(zxw_in)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
c---calculate the couplings as given in Kunszt and Gunion
c---Modified to notation of DKS (ie divided by 2*sw*cw)
c---xw=sin^2 theta_w
      integer:: j
      complex(dp):: zxw_in
      rq1=q1
      rq2=q2
      zxw=zxw_in
      zsin2w=two*sqrt(zxw*(cone-zxw))
      do j=1,nf
      zl(j)=(tau(j)-two*Q(j)*zxw)/zsin2w
      zr(j)=      (-two*Q(j)*zxw)/zsin2w
      enddo

      zle=(-cone-two*(-cone)*zxw)/zsin2w
      zre=(-ctwo*(-cone)*zxw)/zsin2w

      zln=(+cone-two*(+czip)*zxw)/zsin2w
      zrn=czip

      return
      end
