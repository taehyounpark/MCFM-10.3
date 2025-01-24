!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c Matrix element squared for the procoess
c  0 -> H + a(i1) + q(i2) + g(i3) + g(i4)
c with full top mass dependence retained
      use hgggg_integralfill_generic
      use haqgg_pmpp_generic
      use haqgg_pmmp_generic
      use haqgg_pmpm_generic
      use haqgg_pmmm_generic
      use gghgg_dep_params

      implicit none
      include 'Inc/zprods_decl.f'
      include 'Inc/hgggglabels.f'
      include 'Inc/IntResults.f'
      integer i1,i2,i3,i4,h1,h2,h3
      real(dp)::mtsq,mt,mt_yuk,as,vevsq
      real(dp):: msq
      complex(dp) :: amp34(2,2,2),amp43(2,2,2),tmp

      as=real(as_dp,kind=dp)
      vevsq=real(vevsq_dp,kind=dp)
      mt=real(mt_dp,kind=dp)
      mt_yuk=real(mt_yuk_dp,kind=dp)
      mtsq=mt**2

      call hgggg_integralfill(i1,i2,i3,i4,mtsq,Dint,Cint,Bint)

c Ordering of helicities is (quark, gluon, gluon)
c First do (+,-) helicities for (antiquark, quark)
      call haqgg_pmpp(i1,i2,i3,i4,mtsq,za,zb,amp34(1,2,2),Dint,Cint,Bint)
      call haqgg_pmmp(i1,i2,i3,i4,mtsq,za,zb,amp34(1,1,2),Dint,Cint,Bint)
      call haqgg_pmpm(i1,i2,i3,i4,mtsq,za,zb,amp34(1,2,1),Dint,Cint,Bint)
      call haqgg_pmmm(i1,i2,i3,i4,mtsq,za,zb,amp34(1,1,1),Dint,Cint,Bint)

      call haqgg_pmmm(i2,i1,i3,i4,mtsq,zb,za,tmp,Dint,Cint,Bint)
      amp43(1,2,2)=-tmp
      call haqgg_pmpm(i2,i1,i3,i4,mtsq,zb,za,tmp,Dint,Cint,Bint)
      amp43(1,1,2)=-tmp
      call haqgg_pmmp(i2,i1,i3,i4,mtsq,zb,za,tmp,Dint,Cint,Bint)
      amp43(1,2,1)=-tmp
      call haqgg_pmpp(i2,i1,i3,i4,mtsq,zb,za,tmp,Dint,Cint,Bint)
      amp43(1,1,1)=-tmp

c Now call with za<->zb to get (-,+) helicities for (antiquark, quark)
      call haqgg_pmpp(i1,i2,i3,i4,mtsq,zb,za,amp34(2,1,1),Dint,Cint,Bint)
      call haqgg_pmmp(i1,i2,i3,i4,mtsq,zb,za,amp34(2,2,1),Dint,Cint,Bint)
      call haqgg_pmpm(i1,i2,i3,i4,mtsq,zb,za,amp34(2,1,2),Dint,Cint,Bint)
      call haqgg_pmmm(i1,i2,i3,i4,mtsq,zb,za,amp34(2,2,2),Dint,Cint,Bint)

      call haqgg_pmmm(i2,i1,i3,i4,mtsq,za,zb,tmp,Dint,Cint,Bint)
      amp43(2,1,1)=-tmp
      call haqgg_pmpm(i2,i1,i3,i4,mtsq,za,zb,tmp,Dint,Cint,Bint)
      amp43(2,2,1)=-tmp
      call haqgg_pmmp(i2,i1,i3,i4,mtsq,za,zb,tmp,Dint,Cint,Bint)
      amp43(2,1,2)=-tmp
      call haqgg_pmpp(i2,i1,i3,i4,mtsq,za,zb,tmp,Dint,Cint,Bint)
      amp43(2,2,2)=-tmp

      !write(6,*) '|H34 +-++|',abs(amp34(1,2,2))
      !write(6,*) '|H34 +--+|',abs(amp34(1,1,2))
      !write(6,*) '|H34 +-+-|',abs(amp34(1,2,1))

c Add helicities
      msq=zero
      do h1=1,2
      do h2=1,2
      do h3=1,2
      msq=msq+xn*(abs(amp34(h1,h2,h3))**2+abs(amp43(h1,h2,h3))**2)
     &    -one/xn*abs(amp34(h1,h2,h3)+amp43(h1,h2,h3))**2
      enddo
      enddo
      enddo

c Apply overall factor
      msq=msq*V*(as**2*mt*mt_yuk)**2/vevsq

      return


