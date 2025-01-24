!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function AWWjet_qlooptri(j1,j2,j3,j4,j5,j6,j7,za,zb)
c--- This is the amplitude for
c---   0 -> q(p1) q~(p2) nu(p3) e+(p4) mu-(p5) nubar(p6) g(p7)
c--- with helicties q(L), q~(R) and g(R)
c---
c--- It is taken from Eq. (A.14) of arXiv:0710.1832 (CEZ)
c--- and is multiplied by (s127-Mz^2)/s127
c---
c--- Supplemented by Higgs contributions from arXiv:1409.1897, Sec. 3
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
      include 'scalarselect.f'
      include 'yukawas.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      real(dp):: s127,mfsq
      complex(dp):: cF1,cF3,Cmt,B12mt,B1mt,Cmb,B12mb,B1mb,
     & F1mt,F1mb,F3mt,F3mb,prop127,Awwjet_qlooptri,FH

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

c--- Eq. (A.14)
      AWWjet_qlooptri=
     & -cF1*zb(j2,j7)/(4d0*s(j3,j4)*s(j5,j6)*s127)
     & *((za(j1,j3)*zb(j3,j7)+za(j1,j4)*zb(j4,j7)
     &   -za(j1,j5)*zb(j5,j7)-za(j1,j6)*zb(j6,j7))
     &   *za(j3,j5)*zb(j4,j6)
     &  +2d0*(za(j5,j3)*zb(j3,j6)+za(j5,j4)*zb(j4,j6))
     &   *za(j1,j3)*zb(j4,j7)
     &  -2d0*(za(j3,j5)*zb(j5,j4)+za(j3,j6)*zb(j6,j4))
     &   *za(j1,j5)*zb(j6,j7))
     & -cF3*zb(j2,j7)**2/(4d0*s(j1,j2)*s(j3,j4)*s(j5,j6)*s127)
     & *za(j1,j2)*(s(j5,j6)-s(j3,j4))*za(j3,j5)*zb(j4,j6)

c--- apply Z propagator according to scheme
      prop127=cmplx(s127,kind=dp)/cmplx(s127-zmass**2,zmass*zwidth,dp)
      AWWjet_qlooptri=AWWjet_qlooptri*prop127

c--- additional contribution from Higgs diagrams, c.f. 1409.1897, Sec. 3
c Note: we do not account for other finite quark-mass effects in qloop
      FH=2._dp*mt*mt_yuk*(two-(s127-s(j1,j2)-4._dp*mt**2)*Cmt
     &                   +two*s(j1,j2)/(s127-s(j1,j2))*(B12mt-B1mt))
     &  +2._dp*mb*mb_yuk*(two-(s127-s(j1,j2)-4._dp*mb**2)*Cmb
     &                   +two*s(j1,j2)/(s127-s(j1,j2))*(B12mb-B1mb))
      AWWjet_qlooptri=AWWjet_qlooptri
     & -FH/cmplx(s127-hmass**2,hmass*hwidth,dp)/s(j1,j2)
     &  *(za(j1,j2)*zb(j2,j7)**2/(s127-s(j1,j2)))
     &  *za(j3,j5)*zb(j6,j4)/(s(j3,j4)*s(j5,j6))
     & /2._dp

      return
      end

