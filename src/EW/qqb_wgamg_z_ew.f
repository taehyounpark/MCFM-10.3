!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_wgamg_z_ew(p,z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'nwz.f'
      include 'zcouple_cms.f'
      integer:: is
      real(dp):: z,xl12,xl14,xl24,xl16,xl26,xl46,p(mxpart,4),
     & dot,ii_qq,if_qq,fi_qq,ii_qg,ff_qq,tempqg
      real(dp):: QQ1,QQ2

      xl12=log(+two*dot(p,1,2)/musq)
      xl16=log(-two*dot(p,1,6)/musq)
      xl26=log(-two*dot(p,2,6)/musq)
      if (nwz == +1) then
        xl14=log(-two*dot(p,1,4)/musq)
        xl24=log(-two*dot(p,2,4)/musq)
        xl46=log(+two*dot(p,4,6)/musq)
        QQ1=Qu
        QQ2=Qd
      else
        xl14=log(-two*dot(p,1,3)/musq)
        xl24=log(-two*dot(p,2,3)/musq)
        xl46=log(+two*dot(p,3,6)/musq)
        QQ1=Qd
        QQ2=Qu
      endif

c----contributions for one leg
      do is=1,3

c--- These assume we're computing W+ ... need some sort of exchange for W-?
      Q1(q,q,a,is)=abs(zesq)/fourpi/twopi*(
     & QQ1*QQ2*(ii_qq(z,xl12,is))
     &-QQ1*(QQ2-QQ1)*(if_qq(z,xl14,is)+fi_qq(z,xl14,is)))

      Q2(a,a,q,is)=abs(zesq)/fourpi/twopi*(
     & QQ1*QQ2*(ii_qq(z,xl12,is))
     &+QQ2*(QQ2-QQ1)*(if_qq(z,xl24,is)+fi_qq(z,xl24,is)))

      Q1(a,a,q,is)=abs(zesq)/fourpi/twopi*(
     & QQ1*QQ2*(ii_qq(z,xl12,is))
     &+QQ2*(QQ2-QQ1)*(if_qq(z,xl14,is)+fi_qq(z,xl14,is)))

      Q2(q,q,a,is)=abs(zesq)/fourpi/twopi*(
     & QQ1*QQ2*(ii_qq(z,xl12,is))
     &-QQ1*(QQ2-QQ1)*(if_qq(z,xl24,is)+fi_qq(z,xl24,is)))

c--- u + g -> d
      Q1(q,q,g,is)=abs(zesq)/fourpi/twopi*(
     & QQ1*QQ2*(if_qq(z,xl16,is)+fi_qq(z,xl16,is))
     &-QQ1*(QQ2-QQ1)*(if_qq(z,xl14,is)+fi_qq(z,xl14,is))
     &+QQ2*(QQ2-QQ1)*(ff_qq(z,xl46,is)*two))

c--- g + d~ -> u~
      Q2(a,a,g,is)=abs(zesq)/fourpi/twopi*(
     & QQ1*QQ2*(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
     &-QQ1*(QQ2-QQ1)*(ff_qq(z,xl46,is)*two)
     &+QQ2*(QQ2-QQ1)*(if_qq(z,xl24,is)+fi_qq(z,xl24,is)))

c--- d~ + g -> u~
      Q1(a,a,g,is)=abs(zesq)/fourpi/twopi*(
     & QQ1*QQ2*(if_qq(z,xl16,is)+fi_qq(z,xl16,is))
     &-QQ1*(QQ2-QQ1)*(ff_qq(z,xl46,is)*two)
     &+QQ2*(QQ2-QQ1)*(if_qq(z,xl14,is)+fi_qq(z,xl14,is)))

c--- g + u -> d
      Q2(q,q,g,is)=abs(zesq)/fourpi/twopi*(
     & QQ1*QQ2*(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
     &-QQ1*(QQ2-QQ1)*(if_qq(z,xl24,is)+fi_qq(z,xl24,is))
     &+QQ2*(QQ2-QQ1)*(ff_qq(z,xl46,is)*two))



c These are included as part of the Wgaj_a process
c      tempqg=abs(zesq)/fourpi/twopi*xn*ii_qg(z,xl12,is)

c      Q2(a,g,q,is)=tempqg
c      Q2(q,g,a,is)=tempqg
c      Q1(a,g,q,is)=tempqg
c      Q1(q,g,a,is)=tempqg

      enddo

      return
      end
