!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c--- Coefficients of basis integrals in the helicity amplitude for
c----0 -> g+(p1) + g+(p2) + g+(p3) + g+(p4) + h
c--- with full top mass dependence
      use ppppD1x2x34_generic
      use ppppD1x23x4_generic
      use ppppD1x2x3_generic
      use ppppC1x234_generic
      implicit none
      include 'Inc/hgggglabels.f'
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp) :: Dcoeff(dmax),Ccoeff(cmax),Bcoeff(bmax),R(cmax),Rat


      Dcoeff(:)=czip
      Ccoeff(:)=czip
      Bcoeff(:)=czip


      Dcoeff(d1x2x34)=ppppD1x2x34(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x3x41)=ppppD1x2x34(p2,p3,p4,p1,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x4x12)=ppppD1x2x34(p3,p4,p1,p2,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x1x23)=ppppD1x2x34(p4,p1,p2,p3,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d1x4x32)=ppppD1x2x34(p1,p4,p3,p2,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x1x43)=ppppD1x2x34(p2,p1,p4,p3,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x2x14)=ppppD1x2x34(p3,p2,p1,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x3x21)=ppppD1x2x34(p4,p3,p2,p1,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d1x23x4)=ppppD1x23x4(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x34x1)=ppppD1x23x4(p2,p3,p4,p1,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x41x2)=ppppD1x23x4(p3,p4,p1,p2,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x12x3)=ppppD1x23x4(p4,p1,p2,p3,mtsq,za,zb,Cred,I5to4)

      Dcoeff(d1x2x3)=ppppD1x2x3(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d2x3x4)=ppppD1x2x3(p2,p3,p4,p1,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d3x4x1)=ppppD1x2x3(p3,p4,p1,p2,mtsq,za,zb,Cred,I5to4)
      Dcoeff(d4x1x2)=ppppD1x2x3(p4,p1,p2,p3,mtsq,za,zb,Cred,I5to4)

      Ccoeff(c1x234)=ppppC1x234(p1,p2,p3,p4,mtsq,za,zb,R(c1x234))
      Ccoeff(c2x341)=ppppC1x234(p2,p3,p4,p1,mtsq,za,zb,R(c2x341))
      Ccoeff(c3x412)=ppppC1x234(p3,p4,p1,p2,mtsq,za,zb,R(c3x412))
      Ccoeff(c4x123)=ppppC1x234(p4,p1,p2,p3,mtsq,za,zb,R(c4x123))


c      Rat=
c     & (ppppC1x234m2(p1,p2,p3,p4,mtsq,za,zb)
c     & +ppppC1x234m2(p2,p3,p4,p1,mtsq,za,zb)
c     & +ppppC1x234m2(p3,p4,p1,p2,mtsq,za,zb)
c     & +ppppC1x234m2(p4,p1,p2,p3,mtsq,za,zb))/two

      Rat=0.5_dp*(R(c1x234)+R(c2x341)+R(c3x412)+R(c4x123))
      return

