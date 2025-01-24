!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c--- Coefficients of basis integrals in the helicity amplitude for
c----0 -> g+(p1) + g+(p2) + g-(p3) + g-(p4) + h
c--- with full top mass dependence
      use ppmmD1x2x34_generic
      use ppmmD1x4x32_generic
      use ppmmD2x34x1_generic
      use ppmmD1x23x4_generic
      use ppmmD1x2x3_generic
      use ppmmC2x3_generic
      use ppmmC1x23_generic
      use ppmmC23x41_generic
      use ppmmC1x234_generic
      use ppmmB23_generic
      use ppmmB234_generic
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


      Dcoeff(d1x2x34)=ppmmD1x2x34(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x1x43)=ppmmD1x2x34(p2,p1,p4,p3,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x4x12)=ppmmD1x2x34(p3,p4,p1,p2,mtsq,zb,za,Cred,I5to4)
      Dcoeff(d4x3x21)=ppmmD1x2x34(p4,p3,p2,p1,mtsq,zb,za,Cred,I5to4)

      Dcoeff(d1x4x32)=ppmmD1x4x32(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x3x41)=ppmmD1x4x32(p2,p1,p4,p3,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x2x14)=ppmmD1x4x32(p3,p4,p1,p2,mtsq,zb,za,Cred,I5to4)
      Dcoeff(d4x1x23)=ppmmD1x4x32(p4,p3,p2,p1,mtsq,zb,za,Cred,I5to4)

      Dcoeff(d2x34x1)=ppmmD2x34x1(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x12x3)=ppmmD2x34x1(p3,p4,p1,p2,mtsq,zb,za,Cred,I5to4)

      Dcoeff(d1x23x4)=ppmmD1x23x4(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x41x2)=ppmmD1x23x4(p3,p4,p1,p2,mtsq,zb,za,Cred,I5to4)

      Dcoeff(d1x2x3)=ppmmD1x2x3(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x1x2)=ppmmD1x2x3(p2,p1,p4,p3,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x4x1)=ppmmD1x2x3(p3,p4,p1,p2,mtsq,zb,za,Cred,I5to4)
      Dcoeff(d2x3x4)=ppmmD1x2x3(p4,p3,p2,p1,mtsq,zb,za,Cred,I5to4)


      Ccoeff(c2x3)=ppmmC2x3(p1,p2,p3,p4,za,zb)
      Ccoeff(c4x1)=ppmmC2x3(p3,p4,p1,p2,zb,za)

      Ccoeff(c1x23)=ppmmC1x23(p1,p2,p3,p4,za,zb)
      Ccoeff(c2x14)=ppmmC1x23(p2,p1,p4,p3,za,zb)
      Ccoeff(c3x41)=ppmmC1x23(p3,p4,p1,p2,zb,za)
      Ccoeff(c4x32)=ppmmC1x23(p4,p3,p2,p1,zb,za)

      Ccoeff(c23x41)=ppmmC23x41(p1,p2,p3,p4,mtsq,za,zb,R(c23x41))

      Ccoeff(c1x234)=ppmmC1x234(p1,p2,p3,p4,mtsq,za,zb,R(c1x234))
      Ccoeff(c2x341)=ppmmC1x234(p2,p1,p4,p3,mtsq,za,zb,R(c2x341))
      Ccoeff(c3x412)=ppmmC1x234(p3,p4,p1,p2,mtsq,zb,za,R(c3x412))
      Ccoeff(c4x123)=ppmmC1x234(p4,p3,p2,p1,mtsq,zb,za,R(c4x123))


      Bcoeff(b23)=ppmmB23(p1,p2,p3,p4,za,zb)
      Bcoeff(b41)=ppmmB23(p3,p4,p1,p2,zb,za)

      Bcoeff(b234)=ppmmB234(p1,p2,p3,p4,za,zb)
      Bcoeff(b341)=ppmmB234(p2,p1,p4,p3,za,zb)
      Bcoeff(b412)=ppmmB234(p3,p4,p1,p2,zb,za)
      Bcoeff(b123)=ppmmB234(p4,p3,p2,p1,zb,za)

      Bcoeff(b1234)=-Bcoeff(b23)-Bcoeff(b41)
     & -Bcoeff(b234)-Bcoeff(b341)-Bcoeff(b412)-Bcoeff(b123)


c      Rat=
c     & (ppmmC23x41m2(p1,p2,p3,p4,mtsq,za,zb)
c     & +ppmmC1x234m2(p1,p2,p3,p4,mtsq,za,zb)
c     & +ppmmC1x234m2(p2,p1,p4,p3,mtsq,za,zb)
c     & +ppmmC1x234m2(p3,p4,p1,p2,mtsq,zb,za)
c     & +ppmmC1x234m2(p4,p3,p2,p1,mtsq,zb,za))/two

      Rat=0.5_dp*(R(c23x41)+R(c1x234)+R(c2x341)+R(c3x412)+R(c4x123))
      return
