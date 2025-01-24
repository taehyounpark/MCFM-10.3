!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
      function alphas(q,amz,nloop)
          use types
          use LHAPDF
          implicit none
          include 'mxpart.f'
          include 'blha.f'
          real(dp) :: alphas, alphasmz
          real(dp), intent(in) :: q, amz
          integer, intent(in) :: nloop

          if (useblha.ne.0) then
             alphas = alphasmz(q,amz,nloop)
          else
             alphas = getalphas(q)
          endif
      end function
