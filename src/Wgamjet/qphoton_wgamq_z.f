!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qphoton_wgamq_z(p,z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'zcouple_cms.f'
      integer:: is
      real(dp):: z,xl12,p(mxpart,4),dot,ii_qg,tempqg

      xl12=log(+two*dot(p,1,2)/musq)

c----contributions for one leg
      do is=1,3

      tempqg=abs(zesq)/fourpi/twopi*xn*ii_qg(z,xl12,is)

      Q2(a,g,q,is)=tempqg
      Q2(q,g,a,is)=tempqg
      Q1(a,g,q,is)=tempqg
      Q1(q,g,a,is)=tempqg

      enddo

      return
      end
