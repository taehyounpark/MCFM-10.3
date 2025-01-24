!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^+ +H
c     Implementation of Eq.~(4.8) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer::p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp)::ppppD1x2x3_res,ppppE1x2x3x4,ppppE4x1x2x3,fac
      real(dp) s12,s13,s14,s23,s24,s34,s1234
      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      s1234=s12+s13+s14+s23+s24+s34

      fac=(s1234-4._dp*mtsq)/(za(p1,p2)*za(p2,p3)*za(p3,p4)*za(p4,p1))
      ppppE1x2x3x4=mtsq*fac*zb(p1,p2)*za(p2,p3)*zb(p3,p4)*za(p4,p1)
      ppppE4x1x2x3=mtsq*fac*zb(p4,p1)*za(p1,p2)*zb(p2,p3)*za(p3,p4)

      ppppD1x2x3_res=-0.5_dp*fac*s12*s23
     & +Cred(1,I5to4(p4,p1,p2))*ppppE4x1x2x3
     & +Cred(5,I5to4(p1,p2,p3))*ppppE1x2x3x4
      return
