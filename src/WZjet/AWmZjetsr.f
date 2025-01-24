!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWmZjetsr(j1,j2,j3,j4,j5,j6,j7,A7srlo,A7srv)
c--- compute amplitude A^b_7 virtual, including both
c--- leading and subleading pieces, for both positive
c--- and negative gluon helicities, and for both lepton helicities for Z decay
c--- labelling is A7xxx(gluehel, Zlephel)
c--- Note: includes both the usual strong coupling renormalization
c--- and also the finite renormalization in the dred scheme
c--- leading order amplitudes are also returned
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
      include 'nwz.f'
      integer:: j1,j2,j3,j4,j5,j6,j7,iz
      complex(dp):: zlc,zslc,A7sr_slc,A7sr_lc,A7srlo(2,2),A7srv(2,2)
      logical:: swapen

c--- color factors for leading and subleading amplitudes
      zlc=cplx2(3._dp,0._dp)           ! Nc
      zslc=cplx2(-1._dp/3._dp,0._dp)     ! -1/Nc

c--- note: should call AZWWb1 (leading color) second to
c--- get the right sign of Alo for the renormalization below
      verbose=.false.

      if (nwz==-1) then
        swapen=.false.
      else
        swapen=.true.
      endif

      do iz=1,2

      call AWZjetsr_amps('mp',swapen,j1,j2,j3,j4,j5,j6,j7,iz,za,zb,.false.,A7srlo(2,iz),A7sr_lc,A7sr_slc)
      A7srv(2,iz)=A7sr_lc*zlc+A7sr_slc*zslc
      call AWZjetsr_amps('mm',.not.(swapen),j2,j1,j4,j3,j6,j5,j7,iz,zb,za,.true.,A7srlo(1,iz),A7sr_lc,A7sr_slc)
      A7srv(1,iz)=A7sr_lc*zlc+A7sr_slc*zslc

      enddo

c--- strong coupling renormalization in dred scheme
      A7srv(:,:)=A7srv(:,:)-zlc*(b0/xn*epinv-1._dp/6._dp)*A7srlo(:,:)
      if (scheme  ==  'tH-V') then
c--- known translation rules
        A7srv(:,:)=A7srv(:,:)-(Cf+xn/6._dp)*A7srlo(:,:)
      endif

      return
      end

