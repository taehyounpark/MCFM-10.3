!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c--- Coefficients of basis integrals in the helicity amplitude for
c----0 -> g+(p1) + g+(p2) + g+(p3) + g-(p4) + h
c--- with full top mass dependence
      use pppmD1x2x34_generic
      use pppmD1x4x32_generic
      use pppmD2x1x43_generic
      use pppmD2x34x1_generic
      use pppmD4x3x21_generic
      use pppmD1x23x4_generic
      use pppmD2x3x4_generic
      use pppmD1x2x3_generic
      use pppmD3x4x1_generic
      use pppmC3x4_generic
      use pppmC2x34_generic
      use pppmC1x43_generic
      use pppmC4x123_generic
      use pppmC1x234_generic
      use pppmC2x341_generic
      use pppmC12x34_generic
      use pppmB34_generic
      use pppmB234_generic

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


      Dcoeff(d1x2x34)=pppmD1x2x34(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x2x14)=pppmD1x2x34(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d1x4x32)=pppmD1x4x32(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x4x12)=pppmD1x4x32(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d2x1x43)=pppmD2x1x43(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x3x41)=pppmD2x1x43(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d2x34x1)=pppmD2x34x1(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x41x2)=pppmD2x34x1(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d4x3x21)=pppmD4x3x21(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x1x23)=pppmD4x3x21(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d1x23x4)=pppmD1x23x4(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x12x3)=pppmD1x23x4(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d2x3x4)=pppmD2x3x4(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x1x2)=pppmD2x3x4(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d1x2x3)=pppmD1x2x3(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d3x4x1)=pppmD3x4x1(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)


      Ccoeff(c3x4)=pppmC3x4(p1,p2,p3,p4,za,zb)
      Ccoeff(c4x1)=pppmC3x4(p3,p2,p1,p4,za,zb)

      Ccoeff(c2x34)=pppmC2x34(p1,p2,p3,p4,za,zb)
      Ccoeff(c2x14)=pppmC2x34(p3,p2,p1,p4,za,zb)

      Ccoeff(c1x43)=pppmC1x43(p1,p2,p3,p4,za,zb)
      Ccoeff(c3x41)=pppmC1x43(p3,p2,p1,p4,za,zb)

      Ccoeff(c4x123)=pppmC4x123(p1,p2,p3,p4,mtsq,za,zb,R(c4x123))

      Ccoeff(c1x234)=pppmC1x234(p1,p2,p3,p4,mtsq,za,zb,R(c1x234))
      Ccoeff(c3x412)=pppmC1x234(p3,p2,p1,p4,mtsq,za,zb,R(c3x412))

      Ccoeff(c2x341)=pppmC2x341(p1,p2,p3,p4,mtsq,za,zb,R(c2x341))

      Ccoeff(c12x34)=pppmC12x34(p1,p2,p3,p4,mtsq,za,zb,R(c12x34))
      Ccoeff(c23x41)=pppmC12x34(p3,p2,p1,p4,mtsq,za,zb,R(c23x41))


      Bcoeff(b34)=pppmB34(p1,p2,p3,p4,za,zb)
      Bcoeff(b41)=pppmB34(p3,p2,p1,p4,za,zb)

      Bcoeff(b234)=pppmB234(p1,p2,p3,p4,za,zb)
      Bcoeff(b412)=pppmB234(p3,p2,p1,p4,za,zb)
      Bcoeff(b341)=-pppmB234(p2,p3,p1,p4,za,zb)
     & -pppmB234(p2,p1,p3,p4,za,zb)

      Bcoeff(b1234)=-Bcoeff(b34)-Bcoeff(b41)
     & -Bcoeff(b234)-Bcoeff(b412)-Bcoeff(b341)

c      Rat=
c     & (pppmC12x34m2(p1,p2,p3,p4,mtsq,za,zb)
c     & +pppmC12x34m2(p3,p2,p1,p4,mtsq,za,zb)
c     & +pppmC1x234m2(p1,p2,p3,p4,mtsq,za,zb)
c     & +pppmC1x234m2(p3,p2,p1,p4,mtsq,za,zb)
c     & +pppmC4x123m2(p1,p2,p3,p4,mtsq,za,zb)
c     & +pppmC2x341m2(p1,p2,p3,p4,mtsq,za,zb))/two

      Rat=0.5_dp*(R(c12x34)+R(c23x41)
     &           +R(c1x234)+R(c3x412)+R(c4x123)+R(c2x341))
      return


