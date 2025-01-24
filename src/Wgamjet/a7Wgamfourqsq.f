!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a7Wgamfourqsq(j1,j2,j3,j4,j5,j6,j7,za,zb,qqxqqsq)

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      include 'nwz.f'
      include 'qqxqq.f'
c     Performs the square and the sum over polarization of the diagrams
c     which have the W and gamma radiated off an external quark line.
      integer j1,j2,j3,j4,j5,j6,j7,h27,h5,type
      real(dp):: qqxqqsq(12),Qswap(2)
      complex(dp):: upper1(2,2),lower1(2,2),upper2(2,2),lower2(2,2),
     & updown1(2,2),updown2(2,2),downup1(2,2),downup2(2,2),qqxqq(8,2,2)

c     call all basic amplitudes
c     four quark (photon + W on same line (upper1)
c                photon + W on different line (updown1)
      call a7Wgam2q(j1,j2,j3,j4,j5,j6,j7,za,zb,upper1,updown1)
      call a7Wgam2q(j2,j1,j3,j4,j5,j7,j6,za,zb,lower1,downup1)

      if (nwz == -1) then
c     four quark crossed (exhange 7 and 8)
      call a7Wgam2q(j1,j2,j3,j4,j5,j7,j6,za,zb,upper2,updown2)
      call a7Wgam2q(j2,j1,j3,j4,j5,j6,j7,za,zb,lower2,downup2)
      elseif (nwz == +1) then
c    four quark crossed (exchange 1 and 2)
      call a7Wgam2q(j2,j1,j3,j4,j5,j6,j7,za,zb,upper2,updown2)
      call a7Wgam2q(j1,j2,j3,j4,j5,j7,j6,za,zb,lower2,downup2)
      endif

c     qqxqqx(type,h28,h5) etc.
c     as defined in qqxqq.f
c     dc_uc1=1,ds_us1=2,ds_dc1=3,us_uc1=4,
c     dc_uc2=5,ds_us2=6,ds_dc2=7,us_uc2=8,
c     du_uu= 9,dd_ud=10,dd_du=11,ud_uu=12

      if (nwz == -1) then
      Qswap(1)=Q(1)
      Qswap(2)=Q(2)
      elseif (nwz == +1) then
      Qswap(2)=Q(1)
      Qswap(1)=Q(2)
      endif

c     dc_uc (h27:,h5)
      qqxqq(dc_uc1,:,:)=upper1(:,:)+updown1(:,:)*Qswap(2)
c     dcsuc2 (h27,h5) for du_uu (7 and 8 swapped)
      qqxqq(dc_uc2,:,:)=upper2(:,:)+updown2(:,:)*Qswap(2)

c     ds_us  (h27,h5)
      qqxqq(ds_us1,:,:)=upper1(:,:)+updown1(:,:)*Qswap(1)
c     ds_us2 for dd_ud (h16,h5), W on lower line, + 7-8 swap
      qqxqq(ds_us2,:,:)=lower2(:,:)+downup2(:,:)*Qswap(1)

c     ds_dc (h16,h5)
      qqxqq(ds_dc1,:,:)=lower1(:,:)+downup1(:,:)*Qswap(1)
c     ds_dc2 for dd_du (h27,h5),W on upperline, + 7-8 swap
      qqxqq(ds_dc2,:,:)=upper2(:,:)+updown2(:,:)*Qswap(1)

c     us_uc h16,h5
      qqxqq(us_uc1,:,:)=lower1(:,:)+downup1(:,:)*Qswap(2)
c     us_uc2 (h16,h5) for ud_uu W on lower line (7 and 8 swapped)
      qqxqq(us_uc2,:,:)=lower2(:,:)+downup2(:,:)*Qswap(2)

c     square diagrams without identical particle interference
      qqxqqsq(:)=0._dp
      do h5=1,2
      do h27=1,2
      do type=1,8
      qqxqqsq(type)=qqxqqsq(type)+abs(qqxqq(type,h27,h5))**2
      enddo
      enddo
      enddo

c     diagrams with identical particle interference
      do type=9,12
c     Include squared terms from above
      qqxqqsq(type)=qqxqqsq(type-8)+qqxqqsq(type-4)
      do h5=1,2
c     + interference for h27=1;
      qqxqqsq(type)=qqxqqsq(type)
     & +dble(qqxqq(type-8,1,h5)*dconjg(qqxqq(type-4,1,h5)))*2._dp/xn
      enddo
      enddo

      qqxqqsq(:)=aveqq*qqxqqsq(:)

      return
      end
