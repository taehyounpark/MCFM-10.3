!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function getason2pic(musq)
      use LHAPDF, only: getalphas
      use ptveto
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nfl.f'
      real(dp),intent(in) :: musq
      complex(dp) :: getason2pic,a_c
      real(dp) :: as,beta0,beta1
      
      as=getalphas(sqrt(abs(musq)))

      if (timelikemusq) then
! matching scale is time-like, translation to complex alphas from 0809.4283, Eq. (25)
        beta0=11._dp-2._dp*nfl/3._dp
        beta1=102._dp-38*nfl/3._dp
        a_c=beta0*as/4._dp
        getason2pic=(as/(one-im*a_c+beta1/beta0*as/fourpi*log(one-im*a_c)))/twopi
      else
        getason2pic=cmplx(as,0._dp,kind=dp)/twopi
      endif
      
      return
      end

