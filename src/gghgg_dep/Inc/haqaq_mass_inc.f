!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c Matrix element squared for the procoess
c  0 -> H + a(i1) + q(i2) + a(i3) + q(i4)
c with full top mass dependence retained
      use aqaqH_generic
      use gghgg_dep_params
      implicit none
      include 'Inc/zprods_decl.f'
      integer i1,i2,i3,i4,h1,h2
      real(dp)::mtsq,mt,mt_yuk,as,vevsq
      real(dp):: msq,idmsq
      complex(dp) :: amp(2,2),idamp(2,2)

      as=real(as_dp,kind=dp)
      vevsq=real(vevsq_dp,kind=dp)
      mt=real(mt_dp,kind=dp)
      mt_yuk=real(mt_yuk_dp,kind=dp)
      mtsq=mt**2

      call aqaqH(i1,i2,i3,i4,mtsq,za,zb,amp,idamp)
      !write(6,*) '|H4q(+-+-)|',abs(amp(1,1))

c Add helicities
      msq=zero
      idmsq=zero
      do h1=1,2
      do h2=1,2
      msq=msq+abs(amp(h1,h2))**2
      idmsq=idmsq+abs(amp(h1,h2))**2+abs(idamp(h1,h2))**2
      if (h1 == h2) idmsq=idmsq+real(2._dp*amp(h1,h2)*conjg(idamp(h1,h2)))/xn
      enddo
      enddo

c Apply overall factor
      msq=msq*V*(as**2*mt*mt_yuk)**2/vevsq
      idmsq=idmsq*V*(as**2*mt*mt_yuk)**2/vevsq

      return

