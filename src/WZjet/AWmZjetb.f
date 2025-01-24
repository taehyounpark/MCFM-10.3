!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWmZjetb(j1,j2,j3,j4,j5,j6,j7,A7blo,A7bv)
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
      include 'zprods_com.f'
      include 'epinv.f'
      include 'verbose.f'
      include 'b0.f'
      include 'scheme.f'
      integer j1,j2,j3,j4,j5,j6,j7,k5,k6,iz
      complex(dp):: zlc,zslc,A7b_slc,A7b_lc
      complex(dp):: A7blo(2,2),A7bv(2,2)

c--- color factors for leading and subleading amplitudes
      zlc=cplx2(3._dp,0._dp)           ! Nc
      zslc=cplx2(-1._dp/3._dp,0._dp)     ! -1/Nc

c--- note: should call AZWWb1 (leading color) second to
c--- get the right sign of Alo for the renormalization below
      verbose=.false.

      do iz=1,2

      if (iz == 1) then
        k5=j5
        k6=j6
      else
        k5=j6
        k6=j5
      endif

      call AWWjetb_amps('mp',.false.,j1,j2,j3,j4,k5,k6,j7,za,zb,.false.,A7blo(2,iz),A7b_lc,A7b_slc)
      A7bv(2,iz)=A7b_lc*zlc+A7b_slc*zslc
      call AWWjetb_amps('mm',.true.,j2,j1,k5,k6,j3,j4,j7,zb,za,.true.,A7blo(1,iz),A7b_lc,A7b_slc)
      A7bv(1,iz)=A7b_lc*zlc+A7b_slc*zslc

      enddo

c--- strong coupling renormalization in dred scheme
      A7bv(:,:)=A7bv(:,:)-zlc*(b0/xn*epinv-1._dp/6._dp)*A7blo(:,:)
      if (scheme  ==  'tH-V') then
c--- known translation rules
        A7bv(:,:)=A7bv(:,:)-(Cf+xn/6._dp)*A7blo(:,:)
      endif

      return
      end

