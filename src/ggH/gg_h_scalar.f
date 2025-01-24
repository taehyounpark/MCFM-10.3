!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_h_scalar(p,msq)
      use loopI3_generic
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
      include 'hdecaymode.f'
      include 'squark.f'
      include 'scalarselect.f'
      integer:: j,k,iq
      complex(dp)::amp
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      real(dp):: hdecay,gg,Asq,msqhgamgam,mu2,mtsq
      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

c Dummy scale: there is no dependence on this in the result (it is finite)
      mu2=mt**2
c--- set msq=0 to initialize
      msq(:,:)=zip
      amp = czip

      s12=s(1,2)

c   Deal with Higgs decay
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

      do iq=1,3
        mtsq = mloopsq(iq)
        if (iq==1) amp = amp+gloop(iq)*((2*hmass**2-8*mtsq)*loopI3(0d0,0d0,hmass**2,mtsq,mtsq,mtsq,mu2,0)-4d0 )
        if (iq==2 .or. iq==3) amp = amp + gloop(iq)*( -8*mtsq*loopI3(0d0,0d0,hmass**2,mtsq,mtsq,mtsq,mu2,0)-4d0 )
      enddo


      Asq = abs(amp)**2
      gg = as**2/pi**2*Asq

      msq(0,0)=avegg*gg*hdecay
      return
      end
