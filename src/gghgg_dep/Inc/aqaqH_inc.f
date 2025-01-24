!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c     Amplitude for 0->qbar(p1)+q(p2)+Qbar(p3)+Q(p4)+H(p5)
c     Idamp is amplitude for identical quarks
      use aqaqHamp_generic
      implicit none
      include 'Inc/zprods_decl.f'
      integer::p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp)::amp(2,2),idamp(2,2)

      call aqaqHamp(p1,p2,p3,p4,mtsq,za,zb,amp)
      call aqaqHamp(p1,p4,p3,p2,mtsq,za,zb,idamp)

      return


