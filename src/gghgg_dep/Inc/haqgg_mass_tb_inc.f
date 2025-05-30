!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c Matrix element squared for the procoess
c  0 -> H + a(i1) + q(i2) + g(i3) + g(i4)
c with full top, bottom mass dependence retained
c Note that the bottom quark is included with the pole mass in the propagator,
c and a running mass mb(mh) in the Yukawa coupling
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
      integer i1,i2,i3,i4,h1,h2,h3,iq
      real(dp)::mtsq,mt,mt_yuk,mb,mb_yuk,as,vevsq
      real(dp):: msq
      complex(dp) :: amp34(2,2,2),amp43(2,2,2),tmp,
     & amp34t(2,2,2),amp43t(2,2,2),amp34b(2,2,2),amp43b(2,2,2)

      as=real(as_dp,kind=dp)
      vevsq=real(vevsq_dp,kind=dp)
      mt=real(mt_dp,kind=dp)
      mt_yuk=real(mt_yuk_dp,kind=dp)
      mb=real(mb_dp,kind=dp)
      mb_yuk=real(mb_yuk_dp,kind=dp)

      do iq=1,2

      if (iq == 1) mtsq=mt**2
      if (iq == 2) mtsq=mb**2

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

      amp34(:,:,:)=amp34(:,:,:)*mtsq
      amp43(:,:,:)=amp43(:,:,:)*mtsq

      if (iq == 1) then
        amp34t(:,:,:)=amp34(:,:,:)
        amp43t(:,:,:)=amp43(:,:,:)
      endif
      if (iq == 2) then
        amp34b(:,:,:)=amp34(:,:,:)
        amp43b(:,:,:)=amp43(:,:,:)
      endif

      enddo

c      amp34b(:,:,:)=amp34b(:,:,:)*(3.38_dp/mb)       ! Hack for Yukawa with mb(mh)
c      amp43b(:,:,:)=amp43b(:,:,:)*(3.38_dp/mb)       ! Hack for Yukawa with mb(mh)

      amp34t(:,:,:)=amp34t(:,:,:)*(mt_yuk/mt)
      amp43t(:,:,:)=amp43t(:,:,:)*(mt_yuk/mt)
      amp34b(:,:,:)=amp34b(:,:,:)*(mb_yuk/mb)
      amp43b(:,:,:)=amp43b(:,:,:)*(mb_yuk/mb)

      amp34(:,:,:)=amp34t(:,:,:)+amp34b(:,:,:)
      amp43(:,:,:)=amp43t(:,:,:)+amp43b(:,:,:)

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
      msq=msq*V*(as**2)**2/vevsq

      return


