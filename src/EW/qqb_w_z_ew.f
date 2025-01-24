!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_w_z_ew(p,z)
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
      real(dp):: z,xl12,xl14,xl24,p(mxpart,4),dot,ii_qq,if_qq,fi_qq

      xl12=log(+two*dot(p,1,2)/musq)
      xl14=log(-two*dot(p,1,4)/musq)
      xl24=log(-two*dot(p,2,4)/musq)

c----contributions for one leg
      do is=1,3

c--- These assume we're computing W+ ... need some sort of exchange for W-?
      Q1(q,q,a,is)=abs(zesq)/fourpi/twopi*(
     & Qu*Qd*(ii_qq(z,xl12,is))
     &-Qu*Qe*(if_qq(z,xl14,is)+fi_qq(z,xl14,is)))

      Q2(a,a,q,is)=abs(zesq)/fourpi/twopi*(
     & Qu*Qd*(ii_qq(z,xl12,is))
     &+Qd*Qe*(if_qq(z,xl24,is)+fi_qq(z,xl24,is)))

      Q1(a,a,q,is)=abs(zesq)/fourpi/twopi*(
     & Qu*Qd*(ii_qq(z,xl12,is))
     &+Qd*Qe*(if_qq(z,xl14,is)+fi_qq(z,xl14,is)))

      Q2(q,q,a,is)=abs(zesq)/fourpi/twopi*(
     & Qu*Qd*(ii_qq(z,xl12,is))
     &-Qu*Qe*(if_qq(z,xl24,is)+fi_qq(z,xl24,is)))

c These are included as part of the Wgaj_a process
c      tempqg=abs(zesq)/fourpi/twopi*xn*ii_qg(z,xl12,is)

c      Q2(a,g,q,is)=tempqg
c      Q2(q,g,a,is)=tempqg
c      Q1(a,g,q,is)=tempqg
c      Q1(q,g,a,is)=tempqg

      enddo

      return
      end
