!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine gg_h_mass(p,msq)
      implicit none
      include 'types.f'

c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
c----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  b(p3)+b(p4))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      include 'yukawas.f'
      integer:: j,k,n,nmin
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      real(dp):: hdecay,gg,Asq,hprod,msqhgamgam,tn,mq,tn_yuk,mq_yuk
      complex(dp):: hamp,ftn
c      real(dp):: tn,ftn
      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

c--- set msq=0 to initialize
      msq(:,:)=zip

      s12=s(1,2)

      if (hdecaymode == 'none') then
         hdecay = 1._dp
      else
c     Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqhgamgam(s12)
      else
      write(6,*) 'Unimplemented process in gg_hgg_v'
      stop
      endif
      hdecay=hdecay/((s12-hmass**2)**2+(hmass*hwidth)**2)
      endif

      hamp=czip
      nmin=4
      if (mc_yuk == 0 .or. mc == 0) nmin=5
      if (mb_yuk == 0 .or. mb == 0) nmin=6
      do n=nmin,6
      if (n == 4) then
        mq=mc
        mq_yuk=mc_yuk
      endif
      if (n == 5) then
        mq=mb
        mq_yuk=mb_yuk
      endif
      if (n == 6) then
        mq=mt
        mq_yuk=mt_yuk
      endif
      tn = 4._dp*(mq/hmass)**2
      tn_yuk = 4._dp*(mq_yuk/hmass)**2
      if (tn<(1.0_dp)) then
         ftn = 0.5_dp*(log((1._dp+sqrt(1._dp-tn))
     &        /(1._dp-sqrt(1._dp-tn)))-im*pi)**2
      else
         ftn = -2._dp*(asin(1.0_dp/sqrt(tn)))**2
      endif
      hamp=hamp+(3._dp*sqrt(tn*tn_yuk)/4._dp)*(2._dp+(tn-1._dp)*ftn)
      enddo
      hprod=abs(hamp)**2

      Asq=(as/(three*pi))**2/vevsq
      gg=half*Asq*V*s12**2
      msq(0,0)=avegg*gg*hdecay*hprod

      return
      end
