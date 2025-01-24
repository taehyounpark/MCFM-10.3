!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWWjetb(j1,j2,j3,j4,j5,j6,j7,Alomp,Alomm,Alopp,Alopm,A7vmp,A7vmm,A7vpp,A7vpm)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'zprods_com.f'
      include 'verbose.f'
      include 'scheme.f'
      include 'b0.f'
      integer j1,j2,j5,j6,j4,j3,j7
      complex(dp):: zlc,zslc,A7b_slc,A7b_lc,
     & Alomp,Alomm,Alopm,Alopp,A7vmp,A7vmm,A7vpm,A7vpp

c--- color factors for leading and subleading amplitudes
      zlc=cplx2(3._dp,0._dp)           ! Nc
      zslc=cplx2(-1._dp/3._dp,0._dp)     ! -1/Nc

c---- positive-helicity gluon (standard case)

c--- note: should call AZWWb1 (leading color) second to
c--- get the right sign of Alo for the renormalization below
      verbose=.false.
      call AWWjetb_amps('mp',.false.,j1,j2,j3,j4,j5,j6,j7,za,zb,.false.,Alomp,A7b_lc,A7b_slc)
      A7vmp=A7b_lc*zlc+A7b_slc*zslc
      call AWWjetb_amps('mm',.true.,j2,j1,j5,j6,j3,j4,j7,zb,za,.true.,Alomm,A7b_lc,A7b_slc)
      A7vmm=A7b_lc*zlc+A7b_slc*zslc
      call AWWjetb_amps('pp',.true.,j2,j1,j3,j4,j5,j6,j7,za,zb,.false.,Alopp,A7b_lc,A7b_slc)
      A7vpp=A7b_lc*zlc+A7b_slc*zslc
      call AWWjetb_amps('pm',.false.,j1,j2,j5,j6,j3,j4,j7,zb,za,.true.,Alopm,A7b_lc,A7b_slc)
      A7vpm=A7b_lc*zlc+A7b_slc*zslc

c--- strong coupling renormalization in dred scheme
      A7vmp=A7vmp-zlc*(b0/xn*epinv-1._dp/6._dp)*Alomp
      A7vmm=A7vmm-zlc*(b0/xn*epinv-1._dp/6._dp)*Alomm
      A7vpp=A7vpp-zlc*(b0/xn*epinv-1._dp/6._dp)*Alopp
      A7vpm=A7vpm-zlc*(b0/xn*epinv-1._dp/6._dp)*Alopm
      if (scheme  ==  'tH-V') then
c--- known translation rules
        A7vmp=A7vmp-(Cf+xn/6._dp)*Alomp
        A7vmm=A7vmm-(Cf+xn/6._dp)*Alomm
        A7vpp=A7vpp-(Cf+xn/6._dp)*Alopp
        A7vpm=A7vpm-(Cf+xn/6._dp)*Alopm
      endif

      return
      end

