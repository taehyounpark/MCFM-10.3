!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c--- Coefficients of basis integrals in the helicity amplitude for
c----0 -> qb+(p1) + q-(p2) + g-(p3) + g-(p4) + h
c--- with full top mass dependence
      use aqppD4x3x21_generic
      use aqppD3x21x4_generic
      use aqppC12x34_generic
      use aqppC3x12_generic
      use aqppC3x412_generic
      use aqppC4x123_generic
      use aqppB12_generic
      use aqpmppB123_generic
      use aqpmppB412_generic
      implicit none
      include 'Inc/hgggglabels.f'
      include 'Inc/IntResults.f'
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp) :: Dcoeff(dmax),Ccoeff(cmax),Bcoeff(bmax),Rat,Amp,
     & R(cmax)

      Dcoeff(:)=czip
      Ccoeff(:)=czip
      Bcoeff(:)=czip

      Dcoeff(d3x4x12)=-aqppD4x3x21(p2,p1,p4,p3,mtsq,zb,za)
      Dcoeff(d4x3x21)=+aqppD4x3x21(p2,p1,p3,p4,mtsq,zb,za)
      Dcoeff(d4x12x3)=-aqppD3x21x4(p2,p1,p4,p3,mtsq,zb,za)

      Ccoeff(c3x21)=+aqppC3x12(p2,p1,p3,p4,zb,za)
      Ccoeff(c4x12)=-aqppC3x12(p2,p1,p4,p3,zb,za)
      Ccoeff(c12x34)=-aqppC12x34(p2,p1,p4,p3,mtsq,zb,za,R(c12x34))
      Ccoeff(c4x123)=-aqppC3x412(p2,p1,p4,p3,mtsq,zb,za,R(c4x123))
      Ccoeff(c3x412)=-aqppC4x123(p2,p1,p4,p3,mtsq,zb,za,R(c3x412))

      Bcoeff(b12)=aqppB12(p2,p1,p3,p4,zb,za)
      Bcoeff(b412)=-aqpmppB123(p2,p1,p4,p3,zb,za)
      Bcoeff(b123)=-aqpmppB412(p2,p1,p4,p3,zb,za)
      Bcoeff(b1234)=-Bcoeff(b12)-Bcoeff(b123)-Bcoeff(b412)

c      Rat=
c     &-(aqppC12x34m2(p2,p1,p4,p3,mtsq,zb,za)
c     & +aqppC3x412m2(p2,p1,p4,p3,mtsq,zb,za)
c     & +aqppC4x123m2(p2,p1,p4,p3,mtsq,zb,za))/two

      Rat=-0.5_dp*(R(c12x34)+R(c4x123)+R(c3x412))

      Amp=
     & +Dcoeff(d3x4x12)*Dint(6)
     & +Dcoeff(d4x3x21)*Dint(8)
     & +Dcoeff(d4x12x3)*Dint(12)
     & +Ccoeff(c12x34)*Cint(5)
     & +Ccoeff(c3x21)*Cint(11)
     & +Ccoeff(c3x412)*Cint(3)
     & +Ccoeff(c4x12)*Cint(13)
     & +Ccoeff(c3x4)*Cint(17)
     & +Ccoeff(c4x123)*Cint(4)
     & +Bcoeff(b1234)*Bint(9)+Bcoeff(b12)*Bint(5)+Bcoeff(b34)*Bint(7)
     & +Bcoeff(b412)*Bint(4)+Bcoeff(b123)*Bint(1)
     & +Rat

      return

