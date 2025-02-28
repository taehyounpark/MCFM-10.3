!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine breitw(x1,mminsq,mmaxsq,rmass,rwidth,msq,wt)
      implicit none
      include 'types.f'
      include 'constants.f'
c---- Given a number 0<x<1 generate a mass-squared msq and a weight wt
c---- such that mminsq<msq<mmaxsq
c---- points are generated around resonance position rmass, but
c---- breit-wigner should still be included in the matrix element
c     wt is the jacobian between integration in msq and integration in x1
      real(dp):: x1,mminsq,mmaxsq,rmass,rwidth,msq,wt
      real(dp):: almin,almax,al,tanal
      include 'zerowidth.f'

c--- in case the maximum msq is very small, just generate linearly for safety
      if ((mmaxsq < rmass*1.e-3_dp) .and. (zerowidth .eqv. .false.)) then
        msq=mminsq+x1*(mmaxsq-mminsq)
        wt=mmaxsq-mminsq
        return
      endif

      if (zerowidth) then
          tanal=0._dp
          almax=+pi/two
          almin=-pi/two
      else
c--- if generating off-resonance then just do logarithmic sampling in msq
          if (mminsq > rmass**2) then
            wt=one
            call pick(2,msq,mminsq,mmaxsq,x1,wt)
            return
          endif
          almin=atan((mminsq-rmass**2)/rmass/rwidth)
          almax=atan((mmaxsq-rmass**2)/rmass/rwidth)
          al=(almax-almin)*x1+almin
          tanal=tan(al)
      endif

      msq=rmass**2+rmass*rwidth*tanal
c---- bw=(1._dp+tanal**2)*rmass**2*rwidth**2
      wt=(almax-almin)*rmass*rwidth*(1._dp+tanal**2)

      if (msq < 0._dp) then
        msq=mminsq
        wt=0._dp
      endif

      return
      end

