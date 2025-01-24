!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine deltaEW(s,t,ischeme,delta)
c WARNING: due to a bug in QCDLoop2 (qlI3 with complex masses)
c          this routine is incorrect unless scalarselect = 2
c
c     diagrams in units of [a/4/pi]
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'scale.f'
      include 'nf.f'
      include 'zcouple_cms.f'

      real(dp):: s,t,u
      complex(dp):: cw,sw,tot(-2:0),smMWsq,
     & vertexnl(-2:0),vertexdu(-2:0),box(-2:0),WW(-2:0),
     & dZW(-2:0),dZZ(-2:0),dZe(-2:0),dMW2(-2:0),dMZ2(-2:0),dswonsw(-2:0),
     & dZnL(-2:0),dZlL(-2:0),dZuL(-2:0),dZdL(-2:0),dr(-2:0)
      real(dp)::aemon4pi,zmasssq,wmasssq,delta(-2:0)
      real(dp),parameter,dimension(-2:0)::rat=(/0._dp,0._dp,1._dp/)
      real(dp),parameter::a0=1._dp
      real(dp),parameter,dimension(-2:0)::pole=(/0._dp,1._dp,0._dp/)
      integer:: e,ischeme
      complex(dp):: MWsq,MZsq,MHsq,mtsq,swsq,cwsq
      include 'cplx.h'

      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)
      MWsq=cmplx(wmass**2,-wmass*wwidth,kind=dp)
      MHsq=cmplx(hmass**2,-hmass*hwidth,kind=dp)
      mtsq=cmplx(mt**2,-mt*twidth,kind=dp)

      cwsq=MWsq/MZsq
      swsq=1._dp-cwsq

      sw=sqrt(swsq)
      cw=sqrt(cwsq)
      mu=sqrt(musq)
      zmasssq=real(MZsq,kind=dp)
      wmasssq=real(MWsq,kind=dp)
      u=-s-t
      smMWsq=cplx1(s)-MWsq
      aemon4pi=real(zaemmz,dp)/(4._dp*pi)
      call counterterm(ischeme,dZe,dZW,dZZ,dMW2,dMZ2,dswonsw,dZnL,dZlL,d
     & ZuL,dZdL)

!      write(6,*) 'box1',loopI4c(zip,zip,zip,zip,t,s,czip,czip,czip,MWsq,musq,0)
!      write(6,*) 'box2',loopI4c(zip,zip,zip,zip,u,s,czip,czip,czip,MWsq,musq,0)
!      write(6,*) 'box3',loopI4c(zip,zip,zip,zip,t,s,czip,MWsq,czip,MZsq,musq,0)
!      write(6,*) 'box4',loopI4c(zip,zip,zip,zip,u,s,czip,MWsq,czip,MZsq,musq,0)

!      write(6,*) 'tri1 ',loopI3(zip,zip,s,czip,czip,MWsq,musq,0)
!      write(6,*) 'tri2 ',loopI3(zip,zip,s,czip,MZsq,czip,musq,0)
!      write(6,*) 'tri3 ',loopI3(zip,zip,s,MZsq,czip,MWsq,musq,0)
!      write(6,*) 'tri4 ',loopI3(zip,zip,s,czip,czip,czip,musq,0)
!      write(6,*) 'tri5 ',loopI3(zip,zip,s,czip,czip,MWsq,musq,0)
!      write(6,*) 'tri6 ',loopI3(zip,zip,s,czip,MZsq,czip,musq,0)
!      write(6,*) 'tri7 ',loopI3(zip,zip,s,MZsq,czip,MWsq,musq,0)
!      write(6,*) 'tri8 ',loopI3(zip,zip,t,czip,czip,czip,musq,0)
!      write(6,*) 'tri9 ',loopI3(zip,zip,t,czip,MWsq,czip,musq,0)
!      write(6,*) 'tri10',loopI3(zip,zip,t,czip,MZsq,czip,musq,0)
!      write(6,*) 'tri11',loopI3(zip,zip,s,czip,czip,MWsq,musq,0)
!      write(6,*) 'tri12',loopI3(zip,zip,s,MWsq,czip,MZsq,musq,0)

      do e=-2,0
      vertexnl(e)= + loopI2c(zip,czip,czip,musq,e) * ( 2 )
      vertexnl(e) = vertexnl(e) + loopI2c(zip,MWsq,czip,musq,e) * ( 2*
     &    swsq**(-1) + MWsq*swsq**(-1)*s**(-1) )
      vertexnl(e) = vertexnl(e) + loopI2c(zip,MZsq,czip,musq,e) * (  - 2
     &     + 1.D0/2.D0*cwsq**(-2)*MWsq*s**(-1) + cwsq**(-1) - 1.D0/2.D0
     &    *cwsq**(-1)*MWsq*s**(-1) + swsq**(-1) + 1.D0/2.D0*MWsq*
     &    swsq**(-1)*s**(-1) )
      vertexnl(e) = vertexnl(e) + loopI2c(s,czip,czip,musq,e) * (  - 1.D0/
     &    2.D0*cwsq**(-2)*MWsq*s**(-1) - 3.D0/4.D0*cwsq**(-1) + 1.D0/2.D
     &    0*cwsq**(-1)*MWsq*s**(-1) + 3.D0/4.D0*swsq**(-1) + 1.D0/2.D0*
     &    MWsq*swsq**(-1)*s**(-1) )
      vertexnl(e) = vertexnl(e) + loopI2c(s,MWsq,czip,musq,e) * (  - 1 -
     &    MWsq*s**(-1) )
      vertexnl(e) = vertexnl(e) + loopI2c(s,MWsq,MZsq,musq,e) * ( 1 -
     &    swsq**(-1) - 2*MWsq*swsq**(-1)*s**(-1) + MWsq*s**(-1) )
      vertexnl(e) = vertexnl(e) + loopI3(zip,zip,s,czip,czip,MWsq,musq
     & ,e) * ( 2*MWsq )
      vertexnl(e) = vertexnl(e) + loopI3(zip,zip,s,czip,MZsq,czip,musq
     & ,e) * (  - 1.D0/2.D0*cwsq**(-3)*MWsq**2*s**(-1) - cwsq**(-2)*
     &    MWsq + 1.D0/2.D0*cwsq**(-2)*MWsq**2*s**(-1) - 1.D0/2.D0*
     &    cwsq**(-1)*s + cwsq**(-1)*MWsq + 1.D0/2.D0*cwsq**(-1)*MWsq**2
     &    *s**(-1) + 1.D0/2.D0*swsq**(-1)*s + MWsq*swsq**(-1) + 1.D0/2.D
     &    0*MWsq**2*swsq**(-1)*s**(-1) )
      vertexnl(e) = vertexnl(e) + loopI3(zip,zip,s,MZsq,czip,MWsq,musq
     & ,e) * ( 4*MWsq*swsq**(-1) - 2*MWsq + 2*MWsq**2*swsq**(-1)*
     &    s**(-1) )
      vertexnl(e) = vertexnl(e) - 1.D0/2.D0*rat(e)*cwsq**(-1) + 1.D0/2.D
     &    0*rat(e)*swsq**(-1)
      vertexdu(e)= + loopI2c(zip,czip,czip,musq,e) * ( 10.D0/9.D0 )
      vertexdu(e) = vertexdu(e) + loopI2c(zip,MWsq,czip,musq,e) * ( 2*
     &    swsq**(-1) + MWsq*swsq**(-1)*s**(-1) )
      vertexdu(e) = vertexdu(e) + loopI2c(zip,MZsq,czip,musq,e) * (  -
     &    10.D0/9.D0 + 1.D0/18.D0*cwsq**(-2)*MWsq*s**(-1) + 1.D0/9.D0*
     &    cwsq**(-1) - 1.D0/18.D0*cwsq**(-1)*MWsq*s**(-1) + swsq**(-1)
     &     + 1.D0/2.D0*MWsq*swsq**(-1)*s**(-1) )
      vertexdu(e) = vertexdu(e) + loopI2c(s,czip,czip,musq,e) * (  - 1.D0/
     &    18.D0*cwsq**(-2)*MWsq*s**(-1) - 1.D0/12.D0*cwsq**(-1) + 1.D0/
     &    18.D0*cwsq**(-1)*MWsq*s**(-1) + 3.D0/4.D0*swsq**(-1) + 1.D0/2.
     &    D0*MWsq*swsq**(-1)*s**(-1) )
      vertexdu(e) = vertexdu(e) + loopI2c(s,MWsq,czip,musq,e) * (  - 1 -
     &    MWsq*s**(-1) )
      vertexdu(e) = vertexdu(e) + loopI2c(s,MWsq,MZsq,musq,e) * ( 1 -
     &    swsq**(-1) - 2*MWsq*swsq**(-1)*s**(-1) + MWsq*s**(-1) )
      vertexdu(e) = vertexdu(e) + loopI3(zip,zip,s,czip,czip,czip,musq
     & ,e) * ( 4.D0/9.D0*s )
      vertexdu(e) = vertexdu(e) + loopI3(zip,zip,s,czip,czip,MWsq,musq
     & ,e) * ( 2*MWsq )
      vertexdu(e) = vertexdu(e) + loopI3(zip,zip,s,czip,MZsq,czip,musq
     & ,e) * (  - 1.D0/18.D0*cwsq**(-3)*MWsq**2*s**(-1) - 1.D0/9.D0*
     &    cwsq**(-2)*MWsq + 1.D0/18.D0*cwsq**(-2)*MWsq**2*s**(-1) - 1.D0
     &    /18.D0*cwsq**(-1)*s + 1.D0/9.D0*cwsq**(-1)*MWsq + 1.D0/2.D0*
     &    cwsq**(-1)*MWsq**2*s**(-1) + 1.D0/2.D0*swsq**(-1)*s - 4.D0/9.D
     &    0*s + MWsq*swsq**(-1) + 1.D0/2.D0*MWsq**2*swsq**(-1)*s**(-1)
     &     )
      vertexdu(e) = vertexdu(e) + loopI3(zip,zip,s,MZsq,czip,MWsq,musq
     & ,e) * ( 4*MWsq*swsq**(-1) - 2*MWsq + 2*MWsq**2*swsq**(-1)*
     &    s**(-1) )
      vertexdu(e) = vertexdu(e) - 1.D0/18.D0*rat(e)*cwsq**(-1) + 1.D0/2.
     &    D0*rat(e)*swsq**(-1)
      box(e)= + loopI2c(t,czip,czip,musq,e) * (  - 1.D0/3.D0*smMWsq*
     &    cwsq**(-1)*u**(-1) + smMWsq*u**(-1)*swsq**(-1) )
      box(e) = box(e) + loopI2c(s,czip,MWsq,musq,e) * (  - 2.D0/3.D0*
     &    smMWsq*u**(-1) )
      box(e) = box(e) + loopI2c(s,MWsq,MZsq,musq,e) * ( 1.D0/3.D0*smMWsq*
     &    cwsq**(-1)*u**(-1) - smMWsq*u**(-1)*swsq**(-1) + 2.D0/3.D0*
     &    smMWsq*u**(-1) )
      box(e) = box(e) + loopI3(zip,zip,t,czip,czip,czip,musq,e) * (
     &     - 1.D0/3.D0*smMWsq*u**(-2)*s**2 + 1.D0/3.D0*smMWsq*u**(-2)*
     &    MWsq*s - smMWsq*u**(-1)*s + 1.D0/3.D0*smMWsq*u**(-1)*MWsq - 2.
     &    D0/3.D0*smMWsq )
      box(e) = box(e) + loopI3(zip,zip,t,czip,MWsq,czip,musq,e) * (
     &     - 1.D0/6.D0*smMWsq*cwsq**(-2)*u**(-2)*MWsq*s - 1.D0/6.D0*
     &    smMWsq*cwsq**(-2)*u**(-1)*MWsq + 1.D0/6.D0*smMWsq*cwsq**(-1)*
     &    u**(-2)*s**2 + 1.D0/2.D0*smMWsq*cwsq**(-1)*u**(-1)*s + 1.D0/3.
     &    D0*smMWsq*cwsq**(-1) - 1.D0/2.D0*smMWsq*u**(-2)*swsq**(-1)*
     &    s**2 + smMWsq*u**(-2)*MWsq*swsq**(-1)*s - 3.D0/2.D0*smMWsq*
     &    u**(-1)*swsq**(-1)*s + smMWsq*u**(-1)*MWsq*swsq**(-1) -
     &    smMWsq*swsq**(-1) )
      box(e) = box(e) + loopI3(zip,zip,t,czip,MZsq,czip,musq,e) * (
     &     - 1.D0/6.D0*smMWsq*cwsq**(-2)*u**(-2)*MWsq*s - 1.D0/6.D0*
     &    smMWsq*cwsq**(-2)*u**(-1)*MWsq + 1.D0/6.D0*smMWsq*cwsq**(-1)*
     &    u**(-2)*s**2 + 1.D0/2.D0*smMWsq*cwsq**(-1)*u**(-1)*s + 1.D0/3.
     &    D0*smMWsq*cwsq**(-1) - 1.D0/2.D0*smMWsq*u**(-2)*swsq**(-1)*
     &    s**2 + 1.D0/3.D0*smMWsq*u**(-2)*s**2 + smMWsq*u**(-2)*MWsq*
     &    swsq**(-1)*s - 1.D0/3.D0*smMWsq*u**(-2)*MWsq*s - 3.D0/2.D0*
     &    smMWsq*u**(-1)*swsq**(-1)*s + smMWsq*u**(-1)*s + smMWsq*
     &    u**(-1)*MWsq*swsq**(-1) - 1.D0/3.D0*smMWsq*u**(-1)*MWsq -
     &    smMWsq*swsq**(-1) + 2.D0/3.D0*smMWsq )
      box(e) = box(e) + loopI3(zip,zip,s,czip,czip,MWsq,musq,e) * ( 2.D
     &    0/3.D0*smMWsq*u**(-2)*s**2 - 2.D0/3.D0*smMWsq*u**(-2)*MWsq*s
     &     + 4.D0/3.D0*smMWsq*u**(-1)*s + 4*smMWsq )
      box(e) = box(e) + loopI3(zip,zip,s,MWsq,czip,MZsq,musq,e) * ( 1.D
     &    0/3.D0*smMWsq*cwsq**(-2)*u**(-2)*MWsq*s - 1.D0/3.D0*smMWsq*
     &    cwsq**(-1)*u**(-2)*s**2 - 2.D0/3.D0*smMWsq*cwsq**(-1)*u**(-1)
     &    *s + smMWsq*u**(-2)*swsq**(-1)*s**2 - 2.D0/3.D0*smMWsq*
     &    u**(-2)*s**2 - 2*smMWsq*u**(-2)*MWsq*swsq**(-1)*s + 2.D0/3.D0
     &    *smMWsq*u**(-2)*MWsq*s + 2*smMWsq*u**(-1)*swsq**(-1)*s - 4.D0/
     &    3.D0*smMWsq*u**(-1)*s + 4*smMWsq*swsq**(-1) - 4*smMWsq )
      box(e) = box(e) + loopI4c(zip,zip,zip,zip,u,s,czip,czip,czip,
     & MWsq,musq,e) * (  - 4.D0/3.D0*smMWsq*u )
      box(e) = box(e) + loopI4c(zip,zip,zip,zip,u,s,czip,MWsq,czip,
     & MZsq,musq,e) * (  - 1.D0/3.D0*smMWsq*cwsq**(-1)*u - smMWsq*u*
     &    swsq**(-1) + 4.D0/3.D0*smMWsq*u )
      box(e) = box(e) + loopI4c(zip,zip,zip,zip,t,s,czip,czip,czip,
     & MWsq,musq,e) * ( 1.D0/3.D0*smMWsq*u**(-2)*s**3 - 2.D0/3.D0*
     &    smMWsq*u**(-2)*MWsq*s**2 + 1.D0/3.D0*smMWsq*u**(-2)*MWsq**2*s
     &     + smMWsq*u**(-1)*s**2 - 4.D0/3.D0*smMWsq*u**(-1)*MWsq*s + 1.D
     &    0/3.D0*smMWsq*u**(-1)*MWsq**2 + 4.D0/3.D0*smMWsq*s - 2.D0/3.D0
     &    *smMWsq*MWsq + 2.D0/3.D0*smMWsq*u )
      box(e) = box(e) + loopI4c(zip,zip,zip,zip,t,s,czip,MWsq,czip,
     & MZsq,musq,e) * (  - 1.D0/6.D0*smMWsq*cwsq**(-3)*u**(-2)*MWsq**2*
     &    s - 1.D0/6.D0*smMWsq*cwsq**(-3)*u**(-1)*MWsq**2 + 1.D0/3.D0*
     &    smMWsq*cwsq**(-2)*u**(-2)*MWsq*s**2 - 1.D0/6.D0*smMWsq*
     &    cwsq**(-2)*u**(-2)*MWsq**2*s + 2.D0/3.D0*smMWsq*cwsq**(-2)*
     &    u**(-1)*MWsq*s + 1.D0/6.D0*smMWsq*cwsq**(-2)*u**(-1)*MWsq**2
     &     + 1.D0/3.D0*smMWsq*cwsq**(-2)*MWsq - 1.D0/6.D0*smMWsq*
     &    cwsq**(-1)*u**(-2)*s**3 + 2.D0/3.D0*smMWsq*cwsq**(-1)*u**(-2)
     &    *MWsq**2*s - 1.D0/2.D0*smMWsq*cwsq**(-1)*u**(-1)*s**2 + 1.D0/
     &    3.D0*smMWsq*cwsq**(-1)*u**(-1)*MWsq**2 - 2.D0/3.D0*smMWsq*
     &    cwsq**(-1)*s - 1.D0/3.D0*smMWsq*cwsq**(-1)*u + 1.D0/2.D0*
     &    smMWsq*u**(-2)*swsq**(-1)*s**3 - 1.D0/3.D0*smMWsq*u**(-2)*
     &    s**3 - 2*smMWsq*u**(-2)*MWsq*swsq**(-1)*s**2 + 2.D0/3.D0*
     &    smMWsq*u**(-2)*MWsq*s**2 + 2*smMWsq*u**(-2)*MWsq**2*
     &    swsq**(-1)*s - 1.D0/3.D0*smMWsq*u**(-2)*MWsq**2*s + 3.D0/2.D0
     &    *smMWsq*u**(-1)*swsq**(-1)*s**2 )
      box(e) = box(e) + loopI4c(zip,zip,zip,zip,t,s,czip,MWsq,czip,
     & MZsq,musq,e) * (  - smMWsq*u**(-1)*s**2 - 4*smMWsq*u**(-1)*MWsq*
     &    swsq**(-1)*s + 4.D0/3.D0*smMWsq*u**(-1)*MWsq*s + smMWsq*
     &    u**(-1)*MWsq**2*swsq**(-1) - 1.D0/3.D0*smMWsq*u**(-1)*MWsq**2
     &     + 2*smMWsq*swsq**(-1)*s - 4.D0/3.D0*smMWsq*s - 2*smMWsq*MWsq
     &    *swsq**(-1) + 2.D0/3.D0*smMWsq*MWsq + smMWsq*u*swsq**(-1) - 2.
     &    D0/3.D0*smMWsq*u )
      WW(e)= + loopI2c(zip,MHsq,MHsq,musq,e) * ( 1.D0/6.D0*MHsq*
     &    swsq**(-1) )
      WW(e) = WW(e) + loopI2c(zip,MWsq,czip,musq,e) * (  - 2.D0/3.D0*
     &    MWsq**2*s**(-1) )
      WW(e) = WW(e) + loopI2c(zip,MWsq,MHsq,musq,e) * (  - 1.D0/12.D0*
     &    MWsq**2*swsq**(-1)*s**(-1) + 1.D0/6.D0*MHsq*MWsq*swsq**(-1)*
     &    s**(-1) - 1.D0/12.D0*MHsq**2*swsq**(-1)*s**(-1) )
      WW(e) = WW(e) + loopI2c(zip,MWsq,MWsq,musq,e) * ( 5.D0/3.D0*MWsq*
     &    swsq**(-1) )
      WW(e) = WW(e) + loopI2c(zip,MWsq,MZsq,musq,e) * (  - 1.D0/12.D0*
     &    cwsq**(-2)*MWsq**2*s**(-1) - 7.D0/12.D0*cwsq**(-1)*MWsq**2*
     &    s**(-1) + 2.D0/3.D0*MWsq**2*s**(-1) )
      WW(e) = WW(e) + loopI2c(zip,MZsq,MZsq,musq,e) * ( 1.D0/6.D0*
     &    cwsq**(-1)*MWsq + 3.D0/2.D0*MWsq*swsq**(-1) )
      WW(e) = WW(e) + loopI2c(zip,mtsq,czip,musq,e) * ( 1.D0/2.D0*
     &    swsq**(-1)*s**(-1)*mtsq**2 )
      WW(e) = WW(e) + loopI2c(zip,mtsq,mtsq,musq,e) * (  - swsq**(-1)*
     &    mtsq )
      WW(e) = WW(e) + loopI2c(s,czip,czip,musq,e) * ( 3*swsq**(-1)*s )
      WW(e) = WW(e) + loopI2c(s,MWsq,czip,musq,e) * (  - 10.D0/3.D0*s - 4.
     &    D0/3.D0*MWsq + 2.D0/3.D0*MWsq**2*s**(-1) )
      WW(e) = WW(e) + loopI2c(s,MWsq,MHsq,musq,e) * ( 1.D0/12.D0*
     &    swsq**(-1)*s + 5.D0/6.D0*MWsq*swsq**(-1) + 1.D0/12.D0*MWsq**2
     &    *swsq**(-1)*s**(-1) - 1.D0/6.D0*MHsq*swsq**(-1) - 1.D0/6.D0*
     &    MHsq*MWsq*swsq**(-1)*s**(-1) + 1.D0/12.D0*MHsq**2*swsq**(-1)*
     &    s**(-1) )
      WW(e) = WW(e) + loopI2c(s,MWsq,MZsq,musq,e) * ( 1.D0/12.D0*
     &    cwsq**(-2)*MWsq**2*s**(-1) + 5.D0/6.D0*cwsq**(-1)*MWsq + 7.D0/
     &    12.D0*cwsq**(-1)*MWsq**2*s**(-1) - 13.D0/4.D0*swsq**(-1)*s +
     &    10.D0/3.D0*s - 5*MWsq*swsq**(-1) + 4.D0/3.D0*MWsq - 2.D0/3.D0
     &    *MWsq**2*s**(-1) )
      WW(e) = WW(e) + loopI2c(s,mtsq,czip,musq,e) * (  - 1.D0/2.D0*
     &    swsq**(-1)*s**(-1)*mtsq**2 - 1.D0/2.D0*swsq**(-1)*mtsq +
     &    swsq**(-1)*s )
      WW(e) = WW(e) - 13.D0/9.D0*rat(e)*swsq**(-1)*s
      enddo
      vertexnl(:)=aemon4pi*vertexnl(:)
      vertexdu(:)=aemon4pi*vertexdu(:)
      box(:)=aemon4pi*box(:)
      dZZ(:)=aemon4pi*dZZ(:)
      WW(:)=aemon4pi*WW(:)
c The first line of the following expression should also contain -dZW(:)
c and the second line should also contain +dZW(:), hence both are omitted
      tot(:)=
     & -(WW(:)-dMW2(:))/smMWsq
     & +vertexnl(:)+vertexdu(:)
     & +0.5d0*(dZnL(:)+dZlL(:)+dZuL(:)+dZdL(:))
     & +2._dp*(dze(:)-dswonsw(:))
     & +box(:)
      tot(:)=tot(:)+conjg(tot(:)) ! factor of two for square
      delta(:)=real(tot(:),dp)
      return
      end
