!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine fixcms(xx,zz)
      implicit none
      include 'types.f'
c---- fixes coupling in complex mass scheme (zz) according to
c---- its existing real value (xx)
      include 'constants.f'
      include 'nf.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      real(dp):: xx
      complex(dp):: zz
      include 'cplx.h'

      if     (xx == le) then
              zz = zle
      elseif (xx == re) then
              zz = zre
      elseif (xx == ln) then
              zz = zln
      elseif (xx == rn) then
              zz = zrn
      elseif (xx == ln*sqrt(3d0)) then
              zz = zln*sqrt(3d0)
      elseif (xx == rn*sqrt(3d0)) then
              zz = zrn*sqrt(3d0)
      elseif (xx == L(5)*sqrt(3d0)) then
              zz = zL(5)*sqrt(3d0)
      elseif (xx == R(5)*sqrt(3d0)) then
              zz = zR(5)*sqrt(3d0)
      elseif (xx == L(1)*sqrt(3d0*xn)) then
              zz = zL(1)*sqrt(3d0*xn)
      elseif (xx == R(1)*sqrt(3d0*xn)) then
              zz = zR(1)*sqrt(3d0*xn)
      elseif (xx == L(4)*sqrt(2d0*xn)) then
              zz = zL(4)*sqrt(2d0*xn)
      elseif (xx == R(4)*sqrt(2d0*xn)) then
              zz = zR(4)*sqrt(2d0*xn)
      elseif (xx == sqrt(xn*2d0)) then
         zz = cplx1(sqrt(xn*2d0))
      elseif (xx == 1._dp) then
        zz = cone
      elseif (xx == 0._dp) then
        zz = czip
      else
        write(6,*) 'Cannot fix coupling in complex mass scheme: xx = ',xx
        stop
      endif

      return
      end
