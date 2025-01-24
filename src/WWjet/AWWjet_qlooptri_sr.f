!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function AWWjet_qlooptri_sr(j1,j2,j3,j4,j5,j6,j7,za,zb)
c--- This is the amplitude for
c---   0 -> q(p1) q~(p2) nu(p3) e+(p4) mu-(p5) nubar(p6) g(p7)
c--- with helicities q(L), q~(R) and g(R)
c---
c--- where the Z current contains singly-resonant contributions only
c---
c--- It is adapted from Eq. (A.14) of arXiv:0710.1832 (CEZ)
c--- and is multiplied by (s127-Mz^2)/s127
      use loopI2_generic
      use loopI3_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'scale.f'
      include 'nf.f'
      include 'zcouple_cms.f'
      include 'scalarselect.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      real(dp):: s127,mfsq,s34,s56,s345,s346,s356,s456
      complex(dp):: cF1,cF3,Cmt,B12mt,B1mt,Cmb,B12mb,B1mb,
     & F1mt,F1mb,F3mt,F3mb,Awwjet_qlooptri_sr,
     & PW,PZ,MWsq,MZsq,vnL,veL,zab2

c--- statement functions
      zab2(j1,j3,j4,j2)=za(j1,j3)*zb(j3,j2)+za(j1,j4)*zb(j4,j2)
      PW(s56)=s56/(s56-MWsq)
      PZ(s56)=s56/(s56-MZsq)

      s127=s(j1,j2)+s(j1,j7)+s(j2,j7)

c--- relevant scalar integrals
      mfsq=mt**2
      Cmt=loopI3(s(j1,j2),zip,s127,mfsq,mfsq,mfsq,musq,0)
      B12mt=loopI2(s127,mfsq,mfsq,musq,0)
      B1mt=loopI2(s(j1,j2),mfsq,mfsq,musq,0)
      mfsq=mb**2
      Cmb=loopI3(s(j1,j2),zip,s127,mfsq,mfsq,mfsq,musq,0)
      B12mb=loopI2(s127,mfsq,mfsq,musq,0)
      B1mb=loopI2(s(j1,j2),mfsq,mfsq,musq,0)

c--- Eq. (A.7)
      F1mt=1d0/(s(j1,j7)+s(j2,j7))*(2d0+4d0*mt**2*Cmt
     & +(2d0+2d0*s(j1,j2)/(s(j1,j7)+s(j2,j7)))*(B12mt-B1mt))
      F1mb=1d0/(s(j1,j7)+s(j2,j7))*(2d0+4d0*mb**2*Cmb
     & +(2d0+2d0*s(j1,j2)/(s(j1,j7)+s(j2,j7)))*(B12mb-B1mb))

      F3mt=-F1mt+2d0/(s(j1,j7)+s(j2,j7))*(B12mt-B1mt)
      F3mb=-F1mb+2d0/(s(j1,j7)+s(j2,j7))*(B12mb-B1mb)

c--- combination of functions appearing in amplitude
      cF1=F1mt-F1mb
      cF3=F3mt-F3mb

      MWsq=cmplx(wmass**2,-wmass*wwidth,kind=dp)
      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)

      vnL=zln
      veL=zle

      s34=s(j3,j4)
      s56=s(j5,j6)
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s346=s(j3,j4)+s(j3,j6)+s(j4,j6)
      s356=s(j3,j5)+s(j3,j6)+s(j5,j6)
      s456=s(j4,j5)+s(j4,j6)+s(j5,j6)

c overall factor unknown at the moment ...
      AWWjet_qlooptri_sr = PZ(s127)/s127/4._dp * (

     & + PW(s56)/s56*cF3 * (
     &    - za(j3,j5)*zb(j2,j7)**2*zb(j4,j6)/zb(j1,j2)*veL
     &    + za(j3,j5)*zb(j2,j7)**2*zb(j4,j6)/zb(j1,j2)*vnL
     &    )

     & + PW(s56)/s56*cF1 * (
     &    - 2*za(j1,j3)*zb(j2,j7)*zb(j4,j6)*zab2(j5,j6,j4,j7)*vnL/s456
     &    - 2*za(j3,j5)*zb(j2,j7)*zb(j4,j7)*zab2(j1,j5,j3,j6)*veL/s356
     &    )

     & + PW(s34)/s34*cF3 * (
     &    + za(j3,j5)*zb(j2,j7)**2*zb(j4,j6)/zb(j1,j2)*veL
     &    - za(j3,j5)*zb(j2,j7)**2*zb(j4,j6)/zb(j1,j2)*vnL
     &    )

     & + PW(s34)/s34*cF1 * (
     &    + 2*za(j1,j5)*zb(j2,j7)*zb(j4,j6)*zab2(j3,j4,j6,j7)*veL/s346
     &    + 2*za(j3,j5)*zb(j2,j7)*zb(j6,j7)*zab2(j1,j3,j5,j4)*vnL/s345
     &    ) )

c This will get an overall factor of zxw later on, so put in these factors
c to ensure this term comes with the factor sinw/cosw
      AWWjet_qlooptri_sr=AWWjet_qlooptri_sr/sqrt(zxw*(cone-zxw))

      return
      end

c Results for negative helicity gluon below
c      AWWjet_qlooptri_sr =

c     & + PZ(s127)*PW(s34)*s34**-1*cF3*s127**-1 * (
c     &    - 2*za(j1,j7)**2*za(j3,j5)*zb(j4,j6)/za(j1,j2)*veL
c     &    + 2*za(j1,j7)**2*za(j3,j5)*zb(j4,j6)/za(j1,j2)*vnL
c     &    )

c     & + PZ(s127)*PW(s34)*s34**-1*cF1*s127**-1 * (
c     &    - 4*za(j1,j7)*za(j3,j5)*zb(j2,j6)*zab2(j7,j3,j5,j4)*vnL*s345**-1
c     &    - 4*za(j1,j7)*za(j5,j7)*zb(j4,j6)*zab2(j3,j4,j6,j2)*veL*s346**-1
c     &    )

c     & + PZ(s127)*PW(s56)*cF3*s56**-1*s127**-1 * (
c     &    + 2*za(j1,j7)**2*za(j3,j5)*zb(j4,j6)/za(j1,j2)*veL
c     &    - 2*za(j1,j7)**2*za(j3,j5)*zb(j4,j6)/za(j1,j2)*vnL
c     &    )

c     & + PZ(s127)*PW(s56)*cF1*s56**-1*s127**-1 * (
c     &    + 4*za(j1,j7)*za(j3,j5)*zb(j2,j4)*zab2(j7,j5,j3,j6)*veL*s356**-1
c     &    + 4*za(j1,j7)*za(j3,j7)*zb(j4,j6)*zab2(j5,j6,j4,j2)*vnL*s456**-1
c     &    )

