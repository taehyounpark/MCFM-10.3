!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c--- Coefficients of basis integrals in the helicity amplitude for
c----0 -> g+(p1) + g-(p2) + g+(p3) + g-(p4) + h
c--- with full top mass dependence

      use pmpmD4x3x21_generic
      use pmpmD1x23x4_generic
      use pmpmD1x2x3_generic
      use pmpmC3x4_generic
      use pmpmC2x34_generic
      use pmpmC12x34_generic
      use pmpmC1x234_generic
      use pmpmB34_generic
      use pmpmB234_generic
      implicit none
      include 'Inc/hgggglabels.f'
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp) :: Dcoeff(dmax),Ccoeff(cmax),Bcoeff(bmax),Rat,R(cmax)

      Dcoeff(:)=czip
      Ccoeff(:)=czip
      Bcoeff(:)=czip

      Dcoeff(d4x3x21)=pmpmD4x3x21(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x1x23)=pmpmD4x3x21(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x3x41)=pmpmD4x3x21(p1,p4,p3,p2,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x1x43)=pmpmD4x3x21(p3,p4,p1,p2,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d3x4x12)=pmpmD4x3x21(p2,p1,p4,p3,mtsq,zb,za,Cred,I5to4)
      Dcoeff(d3x2x14)=pmpmD4x3x21(p4,p1,p2,p3,mtsq,zb,za,Cred,I5to4)
      Dcoeff(d1x4x32)=pmpmD4x3x21(p2,p3,p4,p1,mtsq,zb,za,Cred,I5to4)
      Dcoeff(d1x2x34)=pmpmD4x3x21(p4,p3,p2,p1,mtsq,zb,za,Cred,I5to4)

      Dcoeff(d1x23x4)=pmpmD1x23x4(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x12x3)=pmpmD1x23x4(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x34x1)=pmpmD1x23x4(p1,p4,p3,p2,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x41x2)=pmpmD1x23x4(p3,p4,p1,p2,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d1x2x3)=pmpmD1x2x3(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x4x1)=pmpmD1x2x3(p3,p4,p1,p2,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x3x4)=pmpmD1x2x3(p2,p3,p4,p1,mtsq,zb,za,Cred,I5to4)
      Dcoeff(d4x1x2)=pmpmD1x2x3(p4,p1,p2,p3,mtsq,zb,za,Cred,I5to4)


      Ccoeff(c3x4)=pmpmC3x4(p1,p2,p3,p4,za,zb)
      Ccoeff(c1x2)=pmpmC3x4(p3,p4,p1,p2,za,zb)
      Ccoeff(c4x1)=pmpmC3x4(p2,p3,p4,p1,zb,za)
      Ccoeff(c2x3)=pmpmC3x4(p4,p1,p2,p3,zb,za)

      Ccoeff(c2x34)=pmpmC2x34(p1,p2,p3,p4,za,zb)
      Ccoeff(c2x14)=pmpmC2x34(p3,p2,p1,p4,za,zb)
      Ccoeff(c4x32)=pmpmC2x34(p1,p4,p3,p2,za,zb)
      Ccoeff(c4x12)=pmpmC2x34(p3,p4,p1,p2,za,zb)

      Ccoeff(c3x41)=pmpmC2x34(p2,p3,p4,p1,zb,za)
      Ccoeff(c3x21)=pmpmC2x34(p4,p3,p2,p1,zb,za)
      Ccoeff(c1x43)=pmpmC2x34(p2,p1,p4,p3,zb,za)
      Ccoeff(c1x23)=pmpmC2x34(p4,p1,p2,p3,zb,za)

      Ccoeff(c12x34)=pmpmC12x34(p1,p2,p3,p4,mtsq,za,zb,R(c12x34))
      Ccoeff(c23x41)=pmpmC12x34(p2,p3,p4,p1,mtsq,zb,za,R(c23x41))

      Ccoeff(c1x234)=pmpmC1x234(p1,p2,p3,p4,mtsq,za,zb,R(c1x234))
      Ccoeff(c3x412)=pmpmC1x234(p3,p4,p1,p2,mtsq,za,zb,R(c3x412))
      Ccoeff(c2x341)=pmpmC1x234(p2,p3,p4,p1,mtsq,zb,za,R(c2x341))
      Ccoeff(c4x123)=pmpmC1x234(p4,p1,p2,p3,mtsq,zb,za,R(c4x123))


      Bcoeff(b34)=pmpmB34(p1,p2,p3,p4,za,zb)
      Bcoeff(b12)=pmpmB34(p3,p4,p1,p2,za,zb)
      Bcoeff(b41)=pmpmB34(p2,p3,p4,p1,zb,za)
      Bcoeff(b23)=pmpmB34(p4,p1,p2,p3,zb,za)

      Bcoeff(b234)=pmpmB234(p1,p2,p3,p4,za,zb)
      Bcoeff(b412)=pmpmB234(p3,p4,p1,p2,za,zb)
      Bcoeff(b341)=pmpmB234(p2,p3,p4,p1,zb,za)
      Bcoeff(b123)=pmpmB234(p4,p1,p2,p3,zb,za)

      Bcoeff(b1234)=-Bcoeff(b34)-Bcoeff(b12)-Bcoeff(b41)-Bcoeff(b23)
     & -Bcoeff(b234)-Bcoeff(b412)-Bcoeff(b341)-Bcoeff(b123)


c      Rat=
c     & (pmpmC12x34m2(p1,p2,p3,p4,mtsq,za,zb)
c     & +pmpmC1x234m2(p1,p2,p3,p4,mtsq,za,zb)
c     & +pmpmC1x234m2(p3,p4,p1,p2,mtsq,za,zb)
c     & +pmpmC12x34m2(p2,p3,p4,p1,mtsq,zb,za)
c     & +pmpmC1x234m2(p4,p1,p2,p3,mtsq,zb,za)
c     & +pmpmC1x234m2(p2,p3,p4,p1,mtsq,zb,za))/two

      Rat=0.5_dp*(R(c12x34)+R(c23x41)+R(c1x234)
     & +R(c3x412)+R(c2x341)+R(c4x123))
      return

