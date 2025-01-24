!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function WWZcurr_sr(he,swapen,j1,jx,j2,jy,za,zb)
      implicit none
      include 'types.f'
c--- this function implements the replacement of a term
c--- of the form <j1 jx> [j2 jy] in the V+1 jet amplitudes
c--- with the expression corresponding to replacing the
c--- current giving |jx> |jy] with the singly-resonant current

c--- note that jx,jy are actually unused here (except for
c--- returning the "usual" value for checking)

c--- the labels for the momenta in the W->WZ process are
c--- given the prefix "p" and are passed in via common block

c Current for W- -> e-(p3) + nu~(p4) + mu-(p5) + mu+(p6)
c muon helicity he
c swapen allows for calls with 3<->4 and 5<->6 to also swap e,n couplings

c note that extra terms are added in order to make this a conserved current, c.f. ZWWcurr_ab.f

      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'masses.f'
      include 'nf.f'
      include 'zcouple_cms.f'
      logical:: swapen
      integer:: j1,jx,j2,jy
      integer:: p1,p2,p3,p4,p5,p6,p7,he
      real(dp)::s34,s56,s3456,s345,s346,s356,s456,qq,qn,qe
      complex(dp)::PW,PZ,MZsq,MWsq,prop56,zab2,vq,veL,vnL,veR
      complex(dp):: WWZcurr_sr
      common/momWWZ/p1,p2,p3,p4,p5,p6,p7
!$omp threadprivate(/momWWZ/)

c--- statement functions
      zab2(j1,jx,jy,j2)=za(j1,jx)*zb(jx,j2)+za(j1,jy)*zb(jy,j2)
      PW(s3456)=s3456/(s3456-MWsq)
      PZ(qq,qe,vq,veL)=qq*qe+vq*veL*prop56
c     end statement functions

      MWsq=cmplx(wmass**2,-wmass*wwidth,kind=dp)
      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)

      qe=-1._dp
      qn=0._dp
      vnL=zln
      veL=zle
      veR=zre

      s34=s(p3,p4)
      s56=s(p5,p6)
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)
      s346=s(p3,p4)+s(p3,p6)+s(p4,p6)
      s356=s(p3,p5)+s(p3,p6)+s(p5,p6)
      s456=s(p4,p5)+s(p4,p6)+s(p5,p6)
      s3456=s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6)
      prop56=s56/(s56-MZsq)

c--- expression we are replacing: za(j1,jx)*zb(j2,jy)

      if (swapen .eqv. .false.) then

        if (he == 1) then
          WWZcurr_sr= PW(s3456)/s56 * (
     &     + PZ(qe,qe,veL,veL) * (
     &       - za(p3,p5)*zb(j2,p4)*zab2(j1,p3,p5,p6)/s356 )
     &     + PZ(qe,qn,veL,vnL) * (
     &       + za(j1,p3)*zb(p4,p6)*zab2(p5,p4,p6,j2)/s456 ) )
c diagrams with two W's
          WWZcurr_sr=WWZcurr_sr + 0.5_dp/zxw*PW(s34)/s34*PW(s3456) * (
     &    - za(j1,p5)*zb(p4,p6)*zab2(p3,p4,p6,j2)/s346 )

          WWZcurr_sr=WWZcurr_sr - 0.5_dp*(
     & + PZ(qe,qe,veL,veL)*PW(s3456)/s56 * (
     &    + 2*za(p3,p5)*zb(p4,p6))
     & + PZ(qe,qn,veL,vnL)*PW(s3456)/s56 * (
     &    - 2*za(p3,p5)*zb(p4,p6))
     & + PW(s34)*PW(s3456) * (
     &    - za(p3,p5)*zb(p4,p6)/s34/zxw)
     &  )*0.5_dp*(zab2(j1,p3,p4,j2)+zab2(j1,p5,p6,j2))
     &  /(s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6))

       else
          WWZcurr_sr= PW(s3456)/s56 * (
     &     + PZ(qe,qe,veR,veL) * (
     &       - za(p3,p6)*zb(j2,p4)*zab2(j1,p3,p6,p5)/s356 )
     &     + PZ(qe,qn,veR,vnL) * (
     &       + za(j1,p3)*zb(p4,p5)*zab2(p6,p4,p5,j2)/s456 ) )

          WWZcurr_sr=WWZcurr_sr - 0.5_dp*(
     &  + PZ(qe,qe,veR,veL)*PW(s3456)/s56 * (
     &     + 2*za(p3,p6)*zb(p4,p5))
     &  + PZ(qe,qn,veR,vnL)*PW(s3456)/s56 * (
     &     - 2*za(p3,p6)*zb(p4,p5))
     &  )*0.5_dp*(zab2(j1,p3,p4,j2)+zab2(j1,p5,p6,j2))
     &  /(s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6))
        endif

      else

        if (he == 1) then
          WWZcurr_sr= PW(s3456)/s56 * (
     &     + PZ(qe,qn,veL,vnL) * (
     &       - za(p3,p5)*zb(j2,p4)*zab2(j1,p3,p5,p6)/s356 )
     &     + PZ(qe,qe,veL,veL) * (
     &       + za(j1,p3)*zb(p4,p6)*zab2(p5,p4,p6,j2)/s456 ) )
c diagrams with two W's
          WWZcurr_sr=WWZcurr_sr - 0.5_dp/zxw*PW(s34)/s34*PW(s3456) * (
     &    - zb(j2,p6)*za(p3,p5)*zab2(j1,p3,p5,p4)/s345 )

          WWZcurr_sr=WWZcurr_sr - 0.5_dp*(
     & + PZ(qe,qn,veL,vnL)*PW(s3456)/s56 * (
     &    + 2*za(p3,p5)*zb(p4,p6))
     & + PZ(qe,qe,veL,veL)*PW(s3456)/s56 * (
     &    - 2*za(p3,p5)*zb(p4,p6))
     & + PW(s34)*PW(s3456) * (
     &    + zb(p4,p6)*za(p3,p5)/s34/zxw)
     &  )*0.5_dp*(zab2(j1,p3,p4,j2)+zab2(j1,p5,p6,j2))
     &  /(s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6))

        else
          WWZcurr_sr= PW(s3456)/s56 * (
     &     + PZ(qe,qn,veR,vnL) * (
     &       - za(p3,p6)*zb(j2,p4)*zab2(j1,p3,p6,p5)/s356 )
     &     + PZ(qe,qe,veR,veL) * (
     &       + za(j1,p3)*zb(p4,p5)*zab2(p6,p4,p5,j2)/s456 ) )

          WWZcurr_sr=WWZcurr_sr - 0.5_dp*(
     &  + PZ(qe,qn,veR,vnL)*PW(s3456)/s56 * (
     &     + 2*za(p3,p6)*zb(p4,p5))
     &  + PZ(qe,qe,veR,veL)*PW(s3456)/s56 * (
     &     - 2*za(p3,p6)*zb(p4,p5))
     &  )*0.5_dp*(zab2(j1,p3,p4,j2)+zab2(j1,p5,p6,j2))
     &  /(s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6))
        endif

      endif

c      ZWWcurr_sr= 2._dp*PW(s56)/s56 * (
c     &  + za(j1,p3)*zb(p4,p6)*zab2(p5,p4,p6,j2)*PZ(qq,qe,vq,veL)/s456
c     &  - za(p3,p5)*zb(j2,p4)*zab2(j1,p3,p5,p6)*PZ(qq,qn,vq,vnL)/s356 )
c     &          + 2._dp*PW(s34)/s34 * (
c     &  + za(j1,p5)*zb(p6,p4)*zab2(p3,p6,p4,j2)*PZ(qq,qn,vq,vnL)/s346
c     &  - za(p5,p3)*zb(j2,p6)*zab2(j1,p5,p3,p4)*PZ(qq,qe,vq,veL)/s345 )

      return
      end

