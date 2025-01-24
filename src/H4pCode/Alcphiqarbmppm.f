!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Alcphiqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Alcphiqarbmppm

      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'deltar.f'
      real(dp):: s12,s13,s41,s23,s24,s34
      real(dp):: s123,s234,s341,s412,mhsq
      complex(dp):: Vlc,Lsm1DS,Lsm1_2me,L0,L1,lnrat,
     & A0phiqarbmppm,A0phidqarbmppm
      integer:: j1,j2,j3,j4

      s12=s(j1,j2)
      s13=s(j1,j3)
      s41=s(j4,j1)
      s23=s(j2,j3)
      s24=s(j2,j4)
      s34=s(j3,j4)
      s123=s12+s13+s23
      s234=s23+s24+s34
      s341=s34+s41+s13
      s412=s41+s12+s24
      mhsq=s12+s13+s41+s23+s24+s34

c---arXIv:09060008v1 Eq.(4.9)
      Vlc = -2._dp*epinv*epinv2
     &  -epinv*lnrat(musq,-s23)-epinv*lnrat(musq,-s41)
     &  -0.5_dp*lnrat(musq,-s23)**2-0.5_dp*lnrat(musq,-s41)**2
     & +13._dp/6._dp*(2._dp*epinv+lnrat(musq,-s12)+lnrat(musq,-s34))
     & - Lsm1_2me(s123,s234,s23,mhsq)
     & - Lsm1_2me(s341,s412,s41,mhsq)
     & - Lsm1DS(s23,s34,s234)
     & - Lsm1DS(s34,s41,s341)
     & - Lsm1DS(s41,s12,s412)
     & - Lsm1DS(s12,s23,s123)
     & + 101._dp/9._dp - deltar/3._dp

c---arXIv:09060008v1 Eq.(4.8)
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      Alcphiqarbmppm=A0phiqarbmppm(j1,j2,j3,j4,za,zb)*Vlc
     & - 0.5_dp*za(j1,j2)*za(j3,j4)*zb(j2,j3)**2
     & *(L1(-s123,-s12)/s12**2+ L1(-s234,-s34)/s34**2)
     & - 2._dp*za(j1,j4)*zb(j2,j3)
     & *(L0(-s123,-s12)/s12+ L0(-s234,-s34)/s34)
     & -2._dp*A0phidqarbmppm(j1,j2,j3,j4,za,zb)

      return
      end
