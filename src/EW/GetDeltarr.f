!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function GetDeltarr()
c Routine not quite working yet ...
c     in units of [alpha/4/pi], but multiplied by [a/2/pi] to give delta
      use loopI2_generic
      use loopI2p_generic
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'zcouple_cms.f'
      include 'scale.f'
      include 'masses.f'
      include 'scalarselect.f'
      real(dp):: GetDeltarr
      complex(dp):: Delr,MWsq,MZsq,MHsq,mtsq,swsq,cwsq

c Match test program
c      wmass=80.357973609878_dp
c      zmass=91.15348061918277_dp
c      mt=173.2_dp
c      hmass=125._dp
c      zaemmz=1d0/137d0
c      musq=zmass**2

      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp) ! zero width for now
      MWsq=cmplx(wmass**2,-wmass*wwidth,kind=dp) ! zero width for now
      MHsq=cmplx(hmass**2,-hmass*hwidth,kind=dp) ! zero width for now
      mtsq=cmplx(mt**2,-mt*twidth,kind=dp) ! zero width for now

      cwsq=MWsq/MZsq
      swsq=1._dp-cwsq

      Delr= ( 26.D0/9.D0 + 6*swsq**(-1) + 7.D0/2.D0*
     &    log(cwsq)*swsq**(-2) - 2*log(cwsq)*swsq**(-1) )
      Delr = Delr + loopI2c(zip,MWsq,czip,musq,0) * (  - 2.D0
     &    /3.D0*swsq**(-1) )
      Delr = Delr + loopI2c(zip,MWsq,MHsq,musq,0) * (  - 1.D0
     &    /12.D0*swsq**(-2) + swsq**(-1) + 1.D0/6.D0*MHsq*MWsq**(-1)*
     &    swsq**(-2) - 1.D0/2.D0*MHsq*MWsq**(-1)*swsq**(-1) - 1.D0/12.D0
     &    *MHsq**2*MWsq**(-2)*swsq**(-2) + 1.D0/6.D0*MHsq**2*MWsq**(-2)
     &    *swsq**(-1) )
      Delr = Delr + loopI2c(zip,MWsq,MWsq,musq,0) * (  - 53.D
     &    0/3.D0 - 4.D0/3.D0*swsq**(-2) + 12*swsq**(-1) + 4*swsq )
      Delr = Delr + loopI2c(zip,MWsq,MZsq,musq,0) * ( 1.D0/
     &    12.D0*cwsq**(-2) + 4.D0/3.D0*cwsq**(-1) - 5*swsq**(-1) )
      Delr = Delr + loopI2c(zip,MZsq,MHsq,musq,0) * ( 1.D0/
     &    12.D0*swsq**(-2) - 1.D0/6.D0*MHsq*MWsq**(-1)*swsq**(-2) + 1.D0
     &    /6.D0*MHsq*MWsq**(-1)*swsq**(-1) + 1.D0/12.D0*MHsq**2*
     &    MWsq**(-2)*swsq**(-2) - 1.D0/6.D0*MHsq**2*MWsq**(-2)*
     &    swsq**(-1) + 1.D0/12.D0*MHsq**2*MWsq**(-2) )
      Delr = Delr + loopI2c(zip,MZsq,MZsq,musq,0) * ( 4.D0/3.
     &    D0*swsq**(-2) - 4.D0/3.D0*swsq**(-1) )
      Delr = Delr + loopI2c(zip,mtsq,czip,musq,0) * ( 1.D0/2.
     &    D0*MWsq**(-2)*swsq**(-2)*mtsq**2 - MWsq**(-2)*swsq**(-1)*
     &    mtsq**2 - 1.D0/2.D0*MWsq**(-1)*swsq**(-1)*mtsq )
      Delr = Delr + loopI2c(zip,mtsq,mtsq,musq,0) * ( 16.D0/
     &    9.D0 - 8.D0/3.D0*MWsq**(-1)*swsq**(-1)*mtsq + 56.D0/9.D0*
     &    MWsq**(-1)*mtsq - 32.D0/9.D0*MWsq**(-1)*swsq*mtsq )
      Delr = Delr + loopI2c(zip,czip,czip,musq,0) * ( 80.D0/9.D
     &    0 )
      Delr = Delr + loopI2c(wmass**2,czip,czip,musq,0) * ( 3*
     &    swsq**(-2) - 6*swsq**(-1) )
      Delr = Delr + loopI2c(wmass**2,MWsq,czip,musq,0) * ( 8 - 4*
     &    swsq**(-1) )
      Delr = Delr + loopI2c(wmass**2,MWsq,MHsq,musq,0) * (
     &    swsq**(-2) - 2*swsq**(-1) - 1.D0/3.D0*MHsq*MWsq**(-1)*
     &    swsq**(-2) + 2.D0/3.D0*MHsq*MWsq**(-1)*swsq**(-1) + 1.D0/12.D0
     &    *MHsq**2*MWsq**(-2)*swsq**(-2) - 1.D0/6.D0*MHsq**2*MWsq**(-2)
     &    *swsq**(-1) )
      Delr = Delr + loopI2c(wmass**2,MWsq,MZsq,musq,0) * (  - 8
     &     - 1.D0/12.D0*cwsq**(-2) - 4.D0/3.D0*cwsq**(-1) - 33.D0/4.D0*
     &    swsq**(-2) + 22*swsq**(-1) )
      Delr = Delr + loopI2c(wmass**2,mtsq,czip,musq,0) * (  - 1.D0
     &    /2.D0*MWsq**(-2)*swsq**(-2)*mtsq**2 + MWsq**(-2)*swsq**(-1)*
     &    mtsq**2 - 1.D0/2.D0*MWsq**(-1)*swsq**(-2)*mtsq + MWsq**(-1)*
     &    swsq**(-1)*mtsq + swsq**(-2) - 2*swsq**(-1) )
      Delr = Delr + loopI2c(zmass**2,czip,czip,musq,0) * (  - 80.D
     &    0/9.D0 - 7.D0/2.D0*swsq**(-2) + 20.D0/3.D0*swsq**(-1) )
      Delr = Delr + loopI2c(zmass**2,MWsq,MWsq,musq,0) * ( 53.D0/
     &    3.D0 + 33.D0/4.D0*swsq**(-2) - 22*swsq**(-1) - 4*swsq )
      Delr = Delr + loopI2c(zmass**2,MZsq,MHsq,musq,0) * (  -
     &    swsq**(-2) + 1.D0/3.D0*MHsq*MWsq**(-1)*swsq**(-2) - 1.D0/3.D0
     &    *MHsq*MWsq**(-1)*swsq**(-1) - 1.D0/12.D0*MHsq**2*MWsq**(-2)*
     &    swsq**(-2) + 1.D0/6.D0*MHsq**2*MWsq**(-2)*swsq**(-1) - 1.D0/
     &    12.D0*MHsq**2*MWsq**(-2) )
      Delr = Delr + loopI2c(zmass**2,mtsq,mtsq,musq,0) * (  - 16.D
     &    0/9.D0 + 1.D0/2.D0*MWsq**(-1)*swsq**(-2)*mtsq + 13.D0/6.D0*
     &    MWsq**(-1)*swsq**(-1)*mtsq - 56.D0/9.D0*MWsq**(-1)*mtsq + 32.D
     &    0/9.D0*MWsq**(-1)*swsq*mtsq - 1.D0/2.D0*swsq**(-2) + 4.D0/3.D0
     &    *swsq**(-1) )

      Delr = Delr + loopI2pc(zip,MWsq,czip,musq,0) * ( 2.D0/3.D0*MWsq )
      Delr = Delr + loopI2pc(zip,MWsq,MHsq,musq,0) * (
     &     - 1.D0/6.D0*MHsq*swsq**(-1) + 1.D0/12.D0*
     &    MHsq**2*MWsq**(-1)*swsq**(-1) + 1.D0/12.D0*MWsq*swsq**(-1)
     &     )
      Delr = Delr + loopI2pc(zip,MWsq,MZsq,musq,0) * (
     &    1.D0/12.D0*MWsq*cwsq**(-2) + 7.D0/12.D0*MWsq*cwsq**(-1)
     &     - 2.D0/3.D0*MWsq )
      Delr = Delr + loopI2pc(zip,mtsq,czip,musq,0) * (  - 1.D
     &    0/2.D0*MWsq**(-1)*swsq**(-1)*mtsq**2 )

      Delr=abs(zesq)/(4._dp*pi)/(4._dp*pi)*Delr

      GetDeltarr=real(Delr,kind=dp)

      return
      end
