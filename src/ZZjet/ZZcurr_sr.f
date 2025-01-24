!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function ZZcurr_sr(Qid,hq,hl1,hl2,j1,jx,j2,jy,za,zb)
      implicit none
      include 'types.f'
c--- this function implements the replacement of a term
c--- of the form <j1 jx> [j2 jy] in the V+1 jet amplitudes
c--- with the expression corresponding to replacing the
c--- current giving |jx> |jy] with the singly-resonant current

c--- note that jx,jy are actually unused here (except for
c--- returning the "usual" value for checking)

c--- the labels for the momenta in the Z->WW process are
c--- given the prefix "p" and are passed in via common block

c quark flavor Qid and helicity hq
c hl1 helicity of 3, hl2 helicity of 5
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'masses.f'
      include 'nf.f'
      include 'zcouple_cms.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      integer j1,jx,j2,jy
      integer p1,p2,p3,p4,p5,p6,p7,Qid,hq,hl1,hl2
      real(dp)::s34,s56,s3456,s345,s346,s356,s456,qq,sign1,sign2
      complex(dp)::PZ,MZsq,zab2,vq,v1,v2
      complex(dp):: ZZcurr_sr
      common/momWWZ/p1,p2,p3,p4,p5,p6,p7
!$omp threadprivate(/momWWZ/)

c--- statement functions
      zab2(j1,jx,jy,j2)=za(j1,jx)*zb(jx,j2)+za(j1,jy)*zb(jy,j2)
      PZ(qq,q1,vq,v1,s56)=qq*q1+vq*v1*s56/(s56-MZsq)
c      PZ(qq,q1,vq,v1,s56)=s56/(s56-MZsq)
c     end statement functions

      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)

      qq=Q(qid)
      if (hq == 1) then
        vq=zL(Qid)
      else
        vq=zR(Qid)
      endif
      if (hl1 == 1) then
        v1=zl1
      else
        v1=zr1
      endif
      if (hl2 == 1) then
        v2=zl2
      else
        v2=zr2
      endif
      s34=s(p3,p4)
      s56=s(p5,p6)
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)
      s346=s(p3,p4)+s(p3,p6)+s(p4,p6)
      s356=s(p3,p5)+s(p3,p6)+s(p5,p6)
      s456=s(p4,p5)+s(p4,p6)+s(p5,p6)
      s3456=s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6)

      if (hl1==2) then
        sign1=-1._dp
      else
        sign1=1._dp
      endif
      if (hl2==2) then
        sign2=-1._dp
      else
        sign2=+1._dp
      endif

c--- expression we are replacing: za(j1,jx)*zb(j2,jy)

      ZZcurr_sr = (
     & + PZ(qq,q1,vq,v1,s3456)*PZ(q1,q2,v1,v2,s56)/s56 * (
     &    + za(j1,p3)*zb(p4,p6)*zab2(p5,p4,p6,j2)/s456
     &    - za(p3,p5)*zb(j2,p4)*zab2(j1,p3,p5,p6)/s356
     &    )*sign1
     & + PZ(qq,q2,vq,v2,s3456)*PZ(q2,q1,v2,v1,s34)/s34 * (
     &    + za(j1,p5)*zb(p6,p4)*zab2(p3,p6,p4,j2)/s346
     &    - za(p5,p3)*zb(j2,p6)*zab2(j1,p5,p3,p4)/s345
     &    )*sign2 )

      return
      end

