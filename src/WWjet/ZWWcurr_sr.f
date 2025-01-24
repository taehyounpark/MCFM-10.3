!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function ZWWcurr_sr(Qid,hq,j1,jx,j2,jy,za,zb)
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
c swapen flips the role of electron and neutrino, for including swapped diagrams
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'masses.f'
      include 'nf.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      integer j1,jx,j2,jy
      integer p1,p2,p3,p4,p5,p6,p7,Qid,hq
      real(dp)::s34,s56,s3456,s345,s346,s356,s456,qq,qn,qe
      complex(dp)::PW,PZ,MZsq,MWsq,prop3456,zab2,vq,veL,vnL
      complex(dp):: ZWWcurr_sr
      common/momWWZ/p1,p2,p3,p4,p5,p6,p7
!$omp threadprivate(/momWWZ/)

c--- statement functions
      zab2(j1,jx,jy,j2)=za(j1,jx)*zb(jx,j2)+za(j1,jy)*zb(jy,j2)
      PW(s56)=s56/(s56-MWsq)
      PZ(qq,qe,vq,veL)=qq*qe+vq*veL*prop3456
c     end statement functions

      MWsq=cmplx(wmass**2,-wmass*wwidth,kind=dp)
      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)

      qe=-1._dp
      qn=0._dp
      vnL=zln
      veL=zle

      qq=Q(qid)
      if (hq == 1) then
        vq=zL(Qid)
      else
        vq=zR(Qid)
      endif
      s34=s(p3,p4)
      s56=s(p5,p6)
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)
      s346=s(p3,p4)+s(p3,p6)+s(p4,p6)
      s356=s(p3,p5)+s(p3,p6)+s(p5,p6)
      s456=s(p4,p5)+s(p4,p6)+s(p5,p6)
      s3456=s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6)
      prop3456=s3456/(s3456-MZsq)

c--- expression we are replacing: za(j1,jx)*zb(j2,jy)

      ZWWcurr_sr= 2._dp*PW(s56)/s56 * (
     &  + za(j1,p3)*zb(p4,p6)*zab2(p5,p4,p6,j2)*PZ(qq,qe,vq,veL)/s456
     &  - za(p3,p5)*zb(j2,p4)*zab2(j1,p3,p5,p6)*PZ(qq,qn,vq,vnL)/s356 )
     &          + 2._dp*PW(s34)/s34 * (
     &  + za(j1,p5)*zb(p6,p4)*zab2(p3,p6,p4,j2)*PZ(qq,qn,vq,vnL)/s346
     &  - za(p5,p3)*zb(j2,p6)*zab2(j1,p5,p3,p4)*PZ(qq,qe,vq,veL)/s345 )

      return
      end

