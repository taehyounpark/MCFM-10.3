!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qphoton_wgamqg_z(p,z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'zcouple_cms.f'
      include 'qcdcouple.f'
      integer:: is
      real(dp):: z,xl12,p(mxpart,4),dot,ii_qg,tempqg,tempqgqed

      xl12=log(+two*dot(p,1,2)/musq)

c----contributions for one leg
      do is=1,3

      tempqg=ason2pi*tr*ii_qg(z,xl12,is)
      tempqgqed=abs(zesq)/fourpi/twopi*xn*ii_qg(z,xl12,is)

      Q1(q,g,g,is)=tempqg
      Q2(q,g,g,is)=tempqg

      Q1(q,g,a,is)=tempqgqed
      Q2(q,g,a,is)=tempqgqed

      enddo

      return
      end
