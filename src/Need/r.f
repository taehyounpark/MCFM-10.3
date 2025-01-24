!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function rpure(pi,pj)
      use types
      implicit none
      include 'userap.f'

      real(dp) :: rpure
      real(dp), intent(in) :: pi(4),pj(4)

      real(dp) :: pti2, ptj2, ei, ej, r1, r2, biti, bitj, delphi, dely
      real(dp), parameter :: tiny = 1.e-9_dp

      pti2=pi(1)**2+pi(2)**2
      ptj2=pj(1)**2+pj(2)**2

      if (userap) then
c--- use rapidities (not pseudorapidities)
        ei=pi(4)
        ej=pj(4)
      else
        ei=sqrt(pti2+pi(3)**2)
        ej=sqrt(ptj2+pj(3)**2)
      endif

      biti=pi(3)/ei
      bitj=pj(3)/ej

      if  ((abs(1d0+biti) < tiny) .or. (abs(1d0-biti) < tiny)
     &     .or.(abs(1d0+bitj) < tiny) .or. (abs(1d0-bitj) < tiny)) then
c-- set to 100 if any of these is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
          dely=100d0
      else
          r1=(1._dp+biti)*(1._dp-bitj)/((1._dp+bitj)*(1._dp-biti))
          dely=0.5_dp*log(r1)
      endif

      r2= (pi(1)*pj(1)+pi(2)*pj(2))/sqrt(pti2*ptj2)
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)

      rpure = sqrt(dely**2+delphi**2)
      end function

      function r(p,i,j)
      use types
      implicit none
      real(dp):: r
c----calculate the jets separation between p(i) and p(j)
      include 'constants.f'
      include 'mxpart.f'
      include 'userap.f'
      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: i,j
      real(dp):: r1,r2,dely,delphi,ei,ej,pti2,ptj2,
     & biti,bitj
      real(dp), parameter:: tiny=1.e-9_dp

      pti2=p(i,1)**2+p(i,2)**2
      ptj2=p(j,1)**2+p(j,2)**2

      if (userap) then
c--- use rapidities (not pseudorapidities)
        ei=p(i,4)
        ej=p(j,4)
      else
        ei=sqrt(pti2+p(i,3)**2)
        ej=sqrt(ptj2+p(j,3)**2)
      endif

c sanity check in case this is called when a momentum is not filled
      if ( (ei < 1.e-13_dp) .or. (ej < 1.e-13_dp)
     & .or.(pti2 < 1.e-26_dp) .or. (ptj2 < 1.e-26_dp)) then
        r=0._dp
        return
      endif

c      r1= (ei+p(i,3))*(ej-p(j,3))/
c     &   ((ej+p(j,3))*(ei-p(i,3)))
      biti=p(i,3)/ei
      bitj=p(j,3)/ej
      if  ((1d0+biti < tiny) .or. (1d0-biti < tiny)
     & .or.(1d0+bitj < tiny) .or. (1d0-bitj < tiny)) then
c-- set to 100 if any of these is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
        dely=100d0
      else
        r1=(one+biti)*(one-bitj)/((one+bitj)*(one-biti))
        dely=0.5_dp*dlog(r1)
      endif

      r2= (p(i,1)*p(j,1)+p(i,2)*p(j,2))/sqrt(pti2*ptj2)
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)

      r=sqrt(dely**2+delphi**2)

      return
      end

      function delphi(pi,pj)
      use ieee_arithmetic
      use types
      implicit none
      include 'tiny.f'

        real(dp) :: delphi
        real(dp), intent(in) :: pi(4), pj(4)

        real(dp) :: pti2, ptj2

        pti2 = pi(1)**2+pi(2)**2
        ptj2 = pj(1)**2+pj(2)**2
        delphi = (pi(1)*pj(1)+pi(2)*pj(2))/sqrt(pti2*ptj2)

        if ((delphi < -1d0) .and. (delphi > -1d0-1d-12)) then
        delphi = -1d0
        endif

        if ((delphi > 1d0) .and. (delphi < 1d0+1d-12)) then
        delphi = 1d0
        endif

        delphi = acos(delphi)

c       if (ieee_is_nan(delphi)) then
c        write (*,*) (pi(1)*pj(1)+pi(2)*pj(2))/sqrt(pti2*ptj2),
c    &      (pi(1)*pj(1)+pi(2)*pj(2)), ptj2,ptj2
c       endif

      end function

