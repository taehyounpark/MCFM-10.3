!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine counterterm(ischeme,dZe,dZW,dZZ,dMW2,dMZ2,dswonsw,
     & dZnL,dZlL,dZuL,dZdL)
c     counterterms in units of [a/4/pi]
      use loopI2_generic
      use loopI2p_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'scale.f'
      include 'nf.f'
      include 'zcouple_cms.f'
      include 'scalarselect.f'
      complex(dp):: cw,sw,B1lam(-2:0),B1W(-2:0),B1Z(-2:0),
     & dZnL(-2:0),dZlL(-2:0),dZuL(-2:0),dZdL(-2:0),
     & dZW(-2:0),dZZ(-2:0),dZe(-2:0),dMW2(-2:0),dMZ2(-2:0),
     & dswonsw(-2:0),Deltasw(-2:0),dr(-2:0),drMZ(-2:0)
      real(dp)::aemon4pi,ReMZsq,REMWsq
      real(dp),parameter,dimension(-2:0)::rat=(/0._dp,0._dp,1._dp/)
      real(dp),parameter::a0=1._dp
      real(dp),parameter,dimension(-2:0)::pole=(/0._dp,1._dp,0._dp/)
      complex(dp)::B00HH(-2:0)=0,B00W0(-2:0)=0,
     & B00WH(-2:0)=0,B00WW(-2:0)=0,B00WZ(-2:0)=0,
     & B00ZZ(-2:0)=0,B00ZH(-2:0)=0,
     & B00t0(-2:0)=0,B00tt(-2:0)=0,
     & B0Z00(-2:0)=0,
     & B0RW00(-2:0)=0,B0RWW0(-2:0)=0,B0RWWH(-2:0)=0,B0RWWZ(-2:0)=0,B0RWt
     & 0(-2:0)=0,
     & B0RZ00(-2:0)=0,B0RZWW(-2:0)=0,B0RZZH(-2:0)=0,B0RZtt(-2:0)=0,
     & B0pRW00(-2:0)=0,B0pRWW0(-2:0)=0,B0pRWWH(-2:0)=0,B0pRWWZ(-2:0)=0,
     & B0pRWt0(-2:0)=0,B0pRZ00(-2:0)=0,B0pRZWW(-2:0)=0,B0pRZZH(-2:0)=0,
     & B0pRZtt(-2:0)=0,B0p0W0(-2:0)=0,B0p0WH(-2:0)=0,B0p0WZ(-2:0)=0,B0p0
     & t0(-2:0)=0
      integer:: e,ischeme
      complex(dp):: MWsq,MZsq,MHsq,mtsq,swsq,cwsq

      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)
      MWsq=cmplx(wmass**2,-wmass*wwidth,kind=dp)
      MHsq=cmplx(hmass**2,-hmass*hwidth,kind=dp)
      mtsq=cmplx(mt**2,-mt*twidth,kind=dp)

      cwsq=MWsq/MZsq
      swsq=1._dp-cwsq

      sw=sqrt(swsq)
      cw=sqrt(cwsq)
      mu=sqrt(musq)
      ReMZsq=real(MZsq,kind=dp)
      ReMWsq=real(MWsq,kind=dp)
      aemon4pi=real(zaemmz,dp)/(4._dp*pi)
      B1lam(:)=czip
      B1W(-2)=czip
      B1W(-1)=1._dp
      B1W(0)=log(musq/MWsq)-0.5_dp
      B1Z(-2)=czip
      B1Z(-1)=1._dp
      B1Z(0)=log(musq/MZsq)-0.5_dp

c      write(6,*) loopI2c(zmass**2,czip,czip,musq,0)
c      write(6,*) loopI2c(zmass**2,MWsq,MWsq,musq,0)
c      write(6,*) loopI2c(zmass**2,MZsq,MHsq,musq,0)
c      write(6,*) loopI2c(zmass**2,cmplx(wmass**2,0.,dp),MWsq,musq,0)
c      write(6,*) loopI2c(zmass**2,cmplx(zmass**2,0.,dp),MHsq,musq,0)
c      write(6,*) loopI2c(zmass**2,cmplx(wmass**2,0.,dp),cmplx(wmass**2,0.,dp),musq,0)
c      write(6,*) loopI2c(zmass**2,cmplx(zmass**2,0.,dp),MHsq,musq,0)
c      write(6,*) loopI2c(zmass**2,mtsq,mtsq,musq,0)
c      write(6,*) loopI2c(zip,MZsq,MHsq,musq,0)
c      write(6,*) loopI2pc(zmass**2,czip,czip,musq,0)
c      write(6,*) loopI2pc(zmass**2,MWsq,MWsq,musq,0)
c      write(6,*) loopI2pc(zmass**2,MZsq,MHsq,musq,0)
c      write(6,*) loopI2pc(zmass**2,mtsq,mtsq,musq,0)
c      write(6,*) zmass**2,MWsq,MWsq,musq
c      write(6,*) zmass**2,MZsq,MHsq,musq
c      stop

      do e=-1,0
      B00HH(e)=loopI2c(zip,MHsq,MHsq,musq,e)
      enddo
      do e=-1,0
      B00W0(e)=loopI2c(zip,MWsq,czip,musq,e)
      enddo
      do e=-1,0
      B00WH(e)=loopI2c(zip,MWsq,MHsq,musq,e)
      enddo
      do e=-1,0
      B00WW(e)=loopI2c(zip,MWsq,MWsq,musq,e)
      enddo
      do e=-1,0
      B00WZ(e)=loopI2c(zip,MWsq,MZsq,musq,e)
      enddo
      do e=-1,0
      B00ZH(e)=loopI2c(zip,MZsq,MHsq,musq,e)
      enddo
      do e=-1,0
      B00ZZ(e)=loopI2c(zip,MZsq,MZsq,musq,e)
      enddo
      do e=-1,0
      B00t0(e)=loopI2c(zip,mtsq,czip,musq,e)
      enddo
      do e=-1,0
      B00tt(e)=loopI2c(zip,mtsq,mtsq,musq,e)
      enddo
      do e=-1,0
      B0Z00(e)=loopI2c(zmass**2,czip,czip,musq,e)
      enddo
      do e=-1,0
      B0RW00(e)=loopI2c(wmass**2,czip,czip,musq,e)
      enddo
      do e=-1,0
      B0RWW0(e)=loopI2c(wmass**2,MWsq,czip,musq,e)
      enddo
      do e=-1,0
      B0RWWH(e)=loopI2c(wmass**2,MWsq,MHsq,musq,e)
      enddo
      do e=-1,0
      B0RWWZ(e)=loopI2c(wmass**2,MWsq,MZsq,musq,e)
      enddo
      do e=-1,0
      B0RWt0(e)=loopI2c(wmass**2,mtsq,czip,musq,e)
      enddo
      do e=-1,0
      B0RZ00(e)=loopI2c(zmass**2,czip,czip,musq,e)
      enddo
      do e=-1,0
      B0RZWW(e)=loopI2c(zmass**2,MWsq,MWsq,musq,e)
      enddo
      do e=-1,0
      B0RZZH(e)=loopI2c(zmass**2,MZsq,MHsq,musq,e)
      enddo
      do e=-1,0
      B0RZtt(e)=loopI2c(zmass**2,mtsq,mtsq,musq,e)
      enddo
      do e=-1,0
      B0pRW00(e)=loopI2pc(wmass**2,czip,czip,musq,e)
      enddo
      do e=-1,0
      B0pRWW0(e)=loopI2pc(wmass**2,MWsq,czip,musq,e)
      enddo
      do e=-1,0
      B0pRWWH(e)=loopI2pc(wmass**2,MWsq,MHsq,musq,e)
      enddo
      do e=-1,0
      B0pRWWZ(e)=loopI2pc(wmass**2,MWsq,MZsq,musq,e)
      enddo
      do e=-1,0
      B0pRWt0(e)=loopI2pc(wmass**2,mtsq,czip,musq,e)
      enddo
      do e=-1,0
      B0pRZ00(e)=loopI2pc(zmass**2,czip,czip,musq,e)
      enddo
      do e=-1,0
      B0pRZWW(e)=loopI2pc(zmass**2,MWsq,MWsq,musq,e)
      enddo
      do e=-1,0
      B0pRZZH(e)=loopI2pc(zmass**2,MZsq,MHsq,musq,e)
      enddo
      do e=-1,0
      B0pRZtt(e)=loopI2pc(zmass**2,mtsq,mtsq,musq,e)
      enddo
      do e=-1,0
      B0p0W0(e)=loopI2pc(zip,MWsq,czip,musq,e)
      enddo
      do e=-1,0
      B0p0WH(e)=loopI2pc(zip,MWsq,MHsq,musq,e)
      enddo
      do e=-1,0
      B0p0WZ(e)=loopI2pc(zip,MWsq,MZsq,musq,e)
      enddo
      do e=-1,0
      B0p0t0(e)=loopI2pc(zip,mtsq,czip,musq,e)
      enddo
      do e=-2,0
      dZe(e)= + B00WW(e) * (  - 7.D0/2.D0 )
      dZe(e) = dZe(e) + B00tt(e) * ( 8.D0/9.D0 )
      dZe(e) = dZe(e) - 1.D0/3.D0*rat(e)
      dZW(e)= + B00W0(e) * (  - 2.D0/3.D0*ReMWsq**(-2)*MWsq**2 )
      dZW(e) = dZW(e) + B00WH(e) * (  - 1.D0/12.D0*ReMWsq**(-2)*MWsq**2
     &    *swsq**(-1) + 1.D0/6.D0*ReMWsq**(-2)*MHsq*MWsq*swsq**(-1) - 1.
     &    D0/12.D0*ReMWsq**(-2)*MHsq**2*swsq**(-1) )
      dZW(e) = dZW(e) + B00WZ(e) * (  - 1.D0/12.D0*cwsq**(-2)*
     &    ReMWsq**(-2)*MWsq**2 - 7.D0/12.D0*cwsq**(-1)*ReMWsq**(-2)*
     &    MWsq**2 + 2.D0/3.D0*ReMWsq**(-2)*MWsq**2 )
      dZW(e) = dZW(e) + B00t0(e) * ( 1.D0/2.D0*ReMWsq**(-2)*swsq**(-1)*
     &    mtsq**2 )
      dZW(e) = dZW(e) + B0RW00(e) * (  - 3*swsq**(-1) )
      dZW(e) = dZW(e) + B0RWW0(e) * ( 10.D0/3.D0 + 2.D0/3.D0*
     &    ReMWsq**(-2)*MWsq**2 )
      dZW(e) = dZW(e) + B0RWWH(e) * ( 1.D0/12.D0*ReMWsq**(-2)*MWsq**2*
     &    swsq**(-1) - 1.D0/6.D0*ReMWsq**(-2)*MHsq*MWsq*swsq**(-1) + 1.D
     &    0/12.D0*ReMWsq**(-2)*MHsq**2*swsq**(-1) - 1.D0/12.D0*
     &    swsq**(-1) )
      dZW(e) = dZW(e) + B0RWWZ(e) * (  - 10.D0/3.D0 + 1.D0/12.D0*
     &    cwsq**(-2)*ReMWsq**(-2)*MWsq**2 + 7.D0/12.D0*cwsq**(-1)*
     &    ReMWsq**(-2)*MWsq**2 - 2.D0/3.D0*ReMWsq**(-2)*MWsq**2 + 13.D0/
     &    4.D0*swsq**(-1) )
      dZW(e) = dZW(e) + B0RWt0(e) * (  - 1.D0/2.D0*ReMWsq**(-2)*
     &    swsq**(-1)*mtsq**2 - swsq**(-1) )
      dZW(e) = dZW(e) + B0pRW00(e) * (  - 3*ReMWsq*swsq**(-1) )
      dZW(e) = dZW(e) + B0pRWW0(e) * (  - 2.D0/3.D0*ReMWsq**(-1)*
     &    MWsq**2 + 4.D0/3.D0*MWsq + 10.D0/3.D0*ReMWsq )
      dZW(e) = dZW(e) + B0pRWWH(e) * (  - 1.D0/12.D0*ReMWsq**(-1)*
     &    MWsq**2*swsq**(-1) + 1.D0/6.D0*ReMWsq**(-1)*MHsq*MWsq*
     &    swsq**(-1) - 1.D0/12.D0*ReMWsq**(-1)*MHsq**2*swsq**(-1) - 5.D0
     &    /6.D0*MWsq*swsq**(-1) + 1.D0/6.D0*MHsq*swsq**(-1) - 1.D0/12.D0
     &    *ReMWsq*swsq**(-1) )
      dZW(e) = dZW(e) + B0pRWWZ(e) * (  - 1.D0/12.D0*cwsq**(-2)*
     &    ReMWsq**(-1)*MWsq**2 - 7.D0/12.D0*cwsq**(-1)*ReMWsq**(-1)*
     &    MWsq**2 - 5.D0/6.D0*cwsq**(-1)*MWsq + 2.D0/3.D0*ReMWsq**(-1)*
     &    MWsq**2 + 5*MWsq*swsq**(-1) - 4.D0/3.D0*MWsq + 13.D0/4.D0*
     &    ReMWsq*swsq**(-1) - 10.D0/3.D0*ReMWsq )
      dZW(e) = dZW(e) + B0pRWt0(e) * ( 1.D0/2.D0*ReMWsq**(-1)*
     &    swsq**(-1)*mtsq**2 + 1.D0/2.D0*swsq**(-1)*mtsq - ReMWsq*
     &    swsq**(-1) )
      dZW(e) = dZW(e) + 13.D0/9.D0*rat(e)*swsq**(-1)
      dZZ(e)= + B00ZH(e) * (  - 1.D0/12.D0*ReMZsq**(-2)*cwsq**(-3)*
     &    MWsq**2 - 1.D0/12.D0*ReMZsq**(-2)*cwsq**(-2)*MWsq**2 + 1.D0/6.
     &    D0*ReMZsq**(-2)*cwsq**(-2)*MHsq*MWsq - 1.D0/12.D0*
     &    ReMZsq**(-2)*cwsq**(-1)*MWsq**2 + 1.D0/6.D0*ReMZsq**(-2)*
     &    cwsq**(-1)*MHsq*MWsq - 1.D0/12.D0*ReMZsq**(-2)*cwsq**(-1)*
     &    MHsq**2 - 1.D0/12.D0*ReMZsq**(-2)*MWsq**2*swsq**(-1) + 1.D0/6.
     &    D0*ReMZsq**(-2)*MHsq*MWsq*swsq**(-1) - 1.D0/12.D0*
     &    ReMZsq**(-2)*MHsq**2*swsq**(-1) )
      dZZ(e) = dZZ(e) + B0RZ00(e) * ( 80.D0/9.D0 - 103.D0/18.D0*
     &    cwsq**(-1) - 7.D0/2.D0*swsq**(-1) )
      dZZ(e) = dZZ(e) + B0RZWW(e) * (  - 3 - 1.D0/12.D0*cwsq**(-1) + 13.
     &    D0/4.D0*swsq**(-1) )
      dZZ(e) = dZZ(e) + B0RZZH(e) * ( 1.D0/12.D0*ReMZsq**(-2)*
     &    cwsq**(-3)*MWsq**2 + 1.D0/12.D0*ReMZsq**(-2)*cwsq**(-2)*
     &    MWsq**2 - 1.D0/6.D0*ReMZsq**(-2)*cwsq**(-2)*MHsq*MWsq + 1.D0/
     &    12.D0*ReMZsq**(-2)*cwsq**(-1)*MWsq**2 - 1.D0/6.D0*
     &    ReMZsq**(-2)*cwsq**(-1)*MHsq*MWsq + 1.D0/12.D0*ReMZsq**(-2)*
     &    cwsq**(-1)*MHsq**2 + 1.D0/12.D0*ReMZsq**(-2)*MWsq**2*
     &    swsq**(-1) - 1.D0/6.D0*ReMZsq**(-2)*MHsq*MWsq*swsq**(-1) + 1.D
     &    0/12.D0*ReMZsq**(-2)*MHsq**2*swsq**(-1) - 1.D0/12.D0*
     &    cwsq**(-1) - 1.D0/12.D0*swsq**(-1) )
      dZZ(e) = dZZ(e) + B0RZtt(e) * ( 16.D0/9.D0 - 17.D0/18.D0*
     &    cwsq**(-1) - 1.D0/2.D0*swsq**(-1) )
      dZZ(e) = dZZ(e) + B0pRZ00(e) * (  - 103.D0/18.D0*ReMZsq*
     &    cwsq**(-1) - 7.D0/2.D0*ReMZsq*swsq**(-1) + 80.D0/9.D0*ReMZsq
     &     )
      dZZ(e) = dZZ(e) + B0pRZWW(e) * (  - 5.D0/3.D0*cwsq**(-1)*MWsq + 5
     &    *MWsq*swsq**(-1) - 4*MWsq - 1.D0/12.D0*ReMZsq*cwsq**(-1) + 13.
     &    D0/4.D0*ReMZsq*swsq**(-1) - 3*ReMZsq )
      dZZ(e) = dZZ(e) + B0pRZZH(e) * (  - 1.D0/12.D0*ReMZsq**(-1)*
     &    cwsq**(-3)*MWsq**2 - 1.D0/12.D0*ReMZsq**(-1)*cwsq**(-2)*
     &    MWsq**2 + 1.D0/6.D0*ReMZsq**(-1)*cwsq**(-2)*MHsq*MWsq - 1.D0/
     &    12.D0*ReMZsq**(-1)*cwsq**(-1)*MWsq**2 + 1.D0/6.D0*
     &    ReMZsq**(-1)*cwsq**(-1)*MHsq*MWsq - 1.D0/12.D0*ReMZsq**(-1)*
     &    cwsq**(-1)*MHsq**2 - 1.D0/12.D0*ReMZsq**(-1)*MWsq**2*
     &    swsq**(-1) + 1.D0/6.D0*ReMZsq**(-1)*MHsq*MWsq*swsq**(-1) - 1.D
     &    0/12.D0*ReMZsq**(-1)*MHsq**2*swsq**(-1) - 5.D0/6.D0*
     &    cwsq**(-2)*MWsq - 5.D0/6.D0*cwsq**(-1)*MWsq + 1.D0/6.D0*
     &    cwsq**(-1)*MHsq - 5.D0/6.D0*MWsq*swsq**(-1) + 1.D0/6.D0*MHsq*
     &    swsq**(-1) - 1.D0/12.D0*ReMZsq*cwsq**(-1) - 1.D0/12.D0*ReMZsq
     &    *swsq**(-1) )
      dZZ(e) = dZZ(e) + B0pRZtt(e) * (  - 7.D0/18.D0*cwsq**(-1)*mtsq +
     &    1.D0/2.D0*swsq**(-1)*mtsq + 32.D0/9.D0*mtsq - 17.D0/18.D0*
     &    ReMZsq*cwsq**(-1) - 1.D0/2.D0*ReMZsq*swsq**(-1) + 16.D0/9.D0*
     &    ReMZsq )
      dZZ(e) = dZZ(e) - 32.D0/9.D0*rat(e) + 19.D0/9.D0*rat(e)*
     &    cwsq**(-1) + 13.D0/9.D0*rat(e)*swsq**(-1)
      dMZ2(e)= + B00HH(e) * ( 1.D0/6.D0*cwsq**(-1)*MHsq + 1.D0/6.D0*
     &    MHsq*swsq**(-1) )
      dMZ2(e) = dMZ2(e) + B00WW(e) * ( 1.D0/3.D0*cwsq**(-1)*MWsq + 3*
     &    MWsq*swsq**(-1) - 4*MWsq )
      dMZ2(e) = dMZ2(e) + B00ZZ(e) * ( 1.D0/6.D0*cwsq**(-2)*MWsq + 1.D0/
     &    6.D0*cwsq**(-1)*MWsq + 1.D0/6.D0*MWsq*swsq**(-1) )
      dMZ2(e) = dMZ2(e) + B00ZH(e) * (  - 1.D0/12.D0*ReMZsq**(-1)*
     &    cwsq**(-3)*MWsq**2 - 1.D0/12.D0*ReMZsq**(-1)*cwsq**(-2)*
     &    MWsq**2 + 1.D0/6.D0*ReMZsq**(-1)*cwsq**(-2)*MHsq*MWsq - 1.D0/
     &    12.D0*ReMZsq**(-1)*cwsq**(-1)*MWsq**2 + 1.D0/6.D0*
     &    ReMZsq**(-1)*cwsq**(-1)*MHsq*MWsq - 1.D0/12.D0*ReMZsq**(-1)*
     &    cwsq**(-1)*MHsq**2 - 1.D0/12.D0*ReMZsq**(-1)*MWsq**2*
     &    swsq**(-1) + 1.D0/6.D0*ReMZsq**(-1)*MHsq*MWsq*swsq**(-1) - 1.D
     &    0/12.D0*ReMZsq**(-1)*MHsq**2*swsq**(-1) )
      dMZ2(e) = dMZ2(e) + B00tt(e) * (  - 17.D0/9.D0*cwsq**(-1)*mtsq -
     &    swsq**(-1)*mtsq + 32.D0/9.D0*mtsq )
      dMZ2(e) = dMZ2(e) + B0RZ00(e) * ( 103.D0/18.D0*ReMZsq*cwsq**(-1)
     &     + 7.D0/2.D0*ReMZsq*swsq**(-1) - 80.D0/9.D0*ReMZsq )
      dMZ2(e) = dMZ2(e) + B0RZWW(e) * ( 5.D0/3.D0*cwsq**(-1)*MWsq - 5*
     &    MWsq*swsq**(-1) + 4*MWsq + 1.D0/12.D0*ReMZsq*cwsq**(-1) - 13.D
     &    0/4.D0*ReMZsq*swsq**(-1) + 3*ReMZsq )
      dMZ2(e) = dMZ2(e) + B0RZZH(e) * ( 1.D0/12.D0*ReMZsq**(-1)*
     &    cwsq**(-3)*MWsq**2 + 1.D0/12.D0*ReMZsq**(-1)*cwsq**(-2)*
     &    MWsq**2 - 1.D0/6.D0*ReMZsq**(-1)*cwsq**(-2)*MHsq*MWsq + 1.D0/
     &    12.D0*ReMZsq**(-1)*cwsq**(-1)*MWsq**2 - 1.D0/6.D0*
     &    ReMZsq**(-1)*cwsq**(-1)*MHsq*MWsq + 1.D0/12.D0*ReMZsq**(-1)*
     &    cwsq**(-1)*MHsq**2 + 1.D0/12.D0*ReMZsq**(-1)*MWsq**2*
     &    swsq**(-1) - 1.D0/6.D0*ReMZsq**(-1)*MHsq*MWsq*swsq**(-1) + 1.D
     &    0/12.D0*ReMZsq**(-1)*MHsq**2*swsq**(-1) + 5.D0/6.D0*
     &    cwsq**(-2)*MWsq + 5.D0/6.D0*cwsq**(-1)*MWsq - 1.D0/6.D0*
     &    cwsq**(-1)*MHsq + 5.D0/6.D0*MWsq*swsq**(-1) - 1.D0/6.D0*MHsq*
     &    swsq**(-1) + 1.D0/12.D0*ReMZsq*cwsq**(-1) + 1.D0/12.D0*ReMZsq
     &    *swsq**(-1) )
      dMZ2(e) = dMZ2(e) + B0RZtt(e) * ( 7.D0/18.D0*cwsq**(-1)*mtsq - 1.D
     &    0/2.D0*swsq**(-1)*mtsq - 32.D0/9.D0*mtsq + 17.D0/18.D0*ReMZsq
     &    *cwsq**(-1) + 1.D0/2.D0*ReMZsq*swsq**(-1) - 16.D0/9.D0*ReMZsq
     &     )
      dMZ2(e) = dMZ2(e) - 19.D0/9.D0*rat(e)*ReMZsq*cwsq**(-1) - 13.D0/9.
     &    D0*rat(e)*ReMZsq*swsq**(-1) + 32.D0/9.D0*rat(e)*ReMZsq
      dMW2(e)= + B00HH(e) * ( 1.D0/6.D0*MHsq*swsq**(-1) )
      dMW2(e) = dMW2(e) + B00W0(e) * (  - 2.D0/3.D0*ReMWsq**(-1)*
     &    MWsq**2 )
      dMW2(e) = dMW2(e) + B00WH(e) * (  - 1.D0/12.D0*ReMWsq**(-1)*
     &    MWsq**2*swsq**(-1) + 1.D0/6.D0*ReMWsq**(-1)*MHsq*MWsq*
     &    swsq**(-1) - 1.D0/12.D0*ReMWsq**(-1)*MHsq**2*swsq**(-1) )
      dMW2(e) = dMW2(e) + B00WW(e) * ( 5.D0/3.D0*MWsq*swsq**(-1) )
      dMW2(e) = dMW2(e) + B00WZ(e) * (  - 1.D0/12.D0*cwsq**(-2)*
     &    ReMWsq**(-1)*MWsq**2 - 7.D0/12.D0*cwsq**(-1)*ReMWsq**(-1)*
     &    MWsq**2 + 2.D0/3.D0*ReMWsq**(-1)*MWsq**2 )
      dMW2(e) = dMW2(e) + B00ZZ(e) * ( 1.D0/6.D0*cwsq**(-1)*MWsq + 3.D0/
     &    2.D0*MWsq*swsq**(-1) )
      dMW2(e) = dMW2(e) + B00t0(e) * ( 1.D0/2.D0*ReMWsq**(-1)*
     &    swsq**(-1)*mtsq**2 )
      dMW2(e) = dMW2(e) + B00tt(e) * (  - swsq**(-1)*mtsq )
      dMW2(e) = dMW2(e) + B0RW00(e) * ( 3*ReMWsq*swsq**(-1) )
      dMW2(e) = dMW2(e) + B0RWW0(e) * ( 2.D0/3.D0*ReMWsq**(-1)*MWsq**2
     &     - 4.D0/3.D0*MWsq - 10.D0/3.D0*ReMWsq )
      dMW2(e) = dMW2(e) + B0RWWH(e) * ( 1.D0/12.D0*ReMWsq**(-1)*MWsq**2
     &    *swsq**(-1) - 1.D0/6.D0*ReMWsq**(-1)*MHsq*MWsq*swsq**(-1) + 1.
     &    D0/12.D0*ReMWsq**(-1)*MHsq**2*swsq**(-1) + 5.D0/6.D0*MWsq*
     &    swsq**(-1) - 1.D0/6.D0*MHsq*swsq**(-1) + 1.D0/12.D0*ReMWsq*
     &    swsq**(-1) )
      dMW2(e) = dMW2(e) + B0RWWZ(e) * ( 1.D0/12.D0*cwsq**(-2)*
     &    ReMWsq**(-1)*MWsq**2 + 7.D0/12.D0*cwsq**(-1)*ReMWsq**(-1)*
     &    MWsq**2 + 5.D0/6.D0*cwsq**(-1)*MWsq - 2.D0/3.D0*ReMWsq**(-1)*
     &    MWsq**2 - 5*MWsq*swsq**(-1) + 4.D0/3.D0*MWsq - 13.D0/4.D0*
     &    ReMWsq*swsq**(-1) + 10.D0/3.D0*ReMWsq )
      dMW2(e) = dMW2(e) + B0RWt0(e) * (  - 1.D0/2.D0*ReMWsq**(-1)*
     &    swsq**(-1)*mtsq**2 - 1.D0/2.D0*swsq**(-1)*mtsq + ReMWsq*
     &    swsq**(-1) )
      dMW2(e) = dMW2(e) - 13.D0/9.D0*rat(e)*ReMWsq*swsq**(-1)
      dZnL(e)= - 1.D0/4.D0*B1Z(e)*cwsq**(-1) - 1.D0/4.D0*B1Z(e)*
     &    swsq**(-1) - 1.D0/2.D0*B1W(e)*swsq**(-1)
      dZlL(e)= + B1Z(e) - 1.D0/4.D0*B1Z(e)*cwsq**(-1) - 1.D0/4.D0*B1Z(e
     &    )*swsq**(-1) - 1.D0/2.D0*B1W(e)*swsq**(-1) - B1lam(e)
      dZuL(e)= + 4.D0/9.D0*B1Z(e) - 1.D0/36.D0*B1Z(e)*cwsq**(-1) - 1.D0/
     &    4.D0*B1Z(e)*swsq**(-1) - 1.D0/2.D0*B1W(e)*swsq**(-1) - 4.D0/9.
     &    D0*B1lam(e)
      dZdL(e)= + 1.D0/9.D0*B1Z(e) - 1.D0/36.D0*B1Z(e)*cwsq**(-1) - 1.D0/
     &    4.D0*B1Z(e)*swsq**(-1) - 1.D0/2.D0*B1W(e)*swsq**(-1) - 1.D0/9.
     &    D0*B1lam(e)
      dr(e)= + B00HH(e) * ( 1.D0/6.D0*MHsq*MWsq**(-1)*swsq**(-1) )
      dr(e) = dr(e) + B00W0(e) * (  - 4.D0/3.D0 )
      dr(e) = dr(e) + B00WH(e) * ( 5.D0/6.D0*swsq**(-1) - 1.D0/6.D0*
     &    MHsq*MWsq**(-1)*swsq**(-1) )
      dr(e) = dr(e) + B00WW(e) * (  - 7 + 17.D0/3.D0*swsq**(-1) )
      dr(e) = dr(e) + B00WZ(e) * ( 4.D0/3.D0 + 5.D0/6.D0*cwsq**(-1) - 5
     &    *swsq**(-1) )
      dr(e) = dr(e) + B00ZZ(e) * ( 1.D0/6.D0*cwsq**(-1) + 3.D0/2.D0*
     &    swsq**(-1) )
      dr(e) = dr(e) + B00t0(e) * (  - 1.D0/2.D0*MWsq**(-1)*swsq**(-1)*
     &    mtsq )
      dr(e) = dr(e) + B00tt(e) * ( 16.D0/9.D0 - MWsq**(-1)*swsq**(-1)*
     &    mtsq )
      dr(e) = dr(e) + B0p0W0(e) * ( 2.D0/3.D0*MWsq )
      dr(e) = dr(e) + B0p0WH(e) * ( 1.D0/12.D0*MWsq*swsq**(-1) - 1.D0/6.
     &    D0*MHsq*swsq**(-1) + 1.D0/12.D0*MHsq**2*MWsq**(-1)*swsq**(-1)
     &     )
      dr(e) = dr(e) + B0p0WZ(e) * ( 1.D0/12.D0*cwsq**(-2)*MWsq + 7.D0/
     &    12.D0*cwsq**(-1)*MWsq - 2.D0/3.D0*MWsq )
      dr(e) = dr(e) + B0p0t0(e) * (  - 1.D0/2.D0*MWsq**(-1)*swsq**(-1)*
     &    mtsq**2 )
      dr(e) = dr(e) - 2.D0/3.D0*rat(e) + 6*rat(e)*swsq**(-1) + 7.D0/2.D0
     &    *rat(e)*log(MWsq*MZsq**(-1))*swsq**(-2) - 2*rat(e)*log(MWsq*
     &    MZsq**(-1))*swsq**(-1)
      drMZ(e)= + B0Z00(e) * (  - 80.D0/9.D0 )
      drMZ(e) = drMZ(e) + 80.D0/27.D0*rat(e)
      enddo
c      write(6,*) 'dZZ(-1)',dZZ(-1)
c      write(6,*) 'dZZ( 0)',dZZ( 0)
c      stop
      dZe(:)=aemon4pi*dZe(:)
      dZW(:)=aemon4pi*dZW(:)
      dZZ(:)=aemon4pi*dZZ(:)
      dMZ2(:)=aemon4pi*dMZ2(:)
      dMW2(:)=aemon4pi*dMW2(:)
      dZnL(:)=aemon4pi*dZnL(:)
      dZlL(:)=aemon4pi*dZlL(:)
      dZuL(:)=aemon4pi*dZuL(:)
      dZdL(:)=aemon4pi*dZdL(:)
      drMZ(:)=aemon4pi*drMZ(:)
      dr(:)=aemon4pi*dr(:)
      dMZ2(:)=dMZ2(:)-(MZsq-ReMZsq)*dZZ(:)
      dMW2(:)=dMW2(:)-(MWsq-ReMWsq)*(dZW(:)+rat(:)*zaemmz/pi)
      dswonsw(:)=-0.5_dp*cwsq/swsq*(dMW2(:)/MWsq-dMZ2(:)/MZsq)
      Deltasw(:)=dswonsw(:)*sqrt(swsq)
      dr(:)=dr(:)-2*dswonsw(:)-dMW2(:)/MWsq
      if (ischeme == 2) dZe(:)=dZe(:)-0.5_dp*drMZ(:)
      if (ischeme == 3) dZe(:)=dZe(:)-0.5_dp*dr(:)

c      write(6,*) 'dZe',dZe(0)
c      write(6,*) 'dZZ',dZZ(0)
c      write(6,*) 'dMZ2',dMZ2(0)
c      write(6,*) 'dZW',dZW(0)
c      write(6,*) 'dMW2',dMW2(0)
c      write(6,*) 'Deltasw',Deltasw(0)
c      write(6,*) 'dswonsw',dswonsw(0)
c      write(6,*) 'dZnL',dZnl(0)
c      write(6,*) 'dZlL',dZlL(0)
c      write(6,*) 'dZul',dZuL(0)
c      write(6,*) 'dZdL',dZdL(0)

      return
      end
