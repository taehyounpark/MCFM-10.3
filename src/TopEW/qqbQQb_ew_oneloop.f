!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqbQQb_ew_oneloop(corr,s,beta,z)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'masses.f'
      include 'epinv.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'anomcoup.f'
      real(dp):: corr(nf),s,t,u,beta,z,alpha,sigma0,gvt,gat,
     &   xI1,xI2,xI3,D6,sw2,cw2,mw,mz,mh,born,f1,f2,
     &   rz,rw,rb,rh,yphi,ys,ini_corr(nf),dFD(nf),T3(nf),gvq(nf),
     &   gaq(nf),qa(5),bxew(nf,-1:0),bxqcd(nf,-1:0),BB(nf,-1:0),fac
      real(dp) :: db0
      integer ep

      ep = 0
c      musq = (two*mt)**2
      alpha = esq/fourpi
c      alpha = 1._dp/126.3_dp
c      as = gsq/fourpi

      mw = wmass
      mz = zmass
      mh = hmass

      sw2 = xw
c on-shell condition cross-check w/ Doreen
c      sw2 = 1._dp - mw**2/mz**2
      cw2 = 1._dp - sw2

      T3 = (/-0.5_dp,0.5_dp,-0.5_dp,0.5_dp,-0.5_dp/)

      gvq(1:5) = (T3(1:5)-two*sw2*Q(1:5))/two/sqrt(sw2*cw2)
      gaq(1:5) = T3(1:5)/two/sqrt(sw2*cw2)
      gw = 0.5_dp/sqrt(two*sw2)

      gvt = gvq(2)
      gat = gaq(2)

c      as = 1._dp
      sigma0 = pi*as**2/8._dp*(Nc**2-1._dp)/Nc**2*beta/s

      t = -s*(1._dp-beta*z)/two+mt**2
      u = -s*(1._dp+beta*z)/two+mt**2

      rz = mz**2/s
      rb = mb**2/s
      rw = mw**2/s
      rh = mh**2/s

      yphi = mb**2/mw**2
      ys = 1._dp - beta**2 + 4._dp*rb

c --- set BB  = 0 for now
c      BB = 0._dp
      fac = pi*alpha*as*(Nc**2-1._dp)*8._dp*s/(s - mz**2)
      BB(:,0) = - fac*(
     &   + (two - beta**2*(1._dp - z**2))*gvq(:)*gvt
     &   + two*beta*z*gaq(:)*gat )

c      write(6,*) 'BB: ', BB(1:2,0)
c      write(6,*) 'born: ', 4._dp*s/(s - mz**2)*(
c     &   + (two - beta**2*(1._dp - z**2))*gvq(1:2)*gvt
c     &   - two*beta*z*gaq(1:2)*gat )
c      stop
c      BB(:,-1) = - two*fac*(-gvq(:)*gvt + 3._dp*beta*z*gaq(:)*gat)

      dFD = - alpha/8._dp/pi*((gvq**2+gaq**2)*f1(rz)+two*gw**2*f1(rw))

      born = sigma0*(two - beta**2 + beta**2*z**2)

      ini_corr = two*born*dFD

      qa(1) =
     & (0.25_dp*alpha*sigma0*(-two*(gat**2 + gvt**2)*(1._dp + z**2)
     & + two*((1._dp - beta**2)*(-3._dp*gat**2 + gvt**2)
     & + two*(gat**2 + gvt**2)*rz)*s*(two - beta**2*(1._dp
     & - z**2))*db0(mt**2,mz**2,mt**2) + (8._dp*(gat**2
     & + gvt**2)*(1._dp + z**2)*(-xI1(mt**2,musq,ep)
     & + xI1(mz**2,musq,ep)))/((1._dp - beta**2)*s)
     & - 4._dp*(gat**2*(1._dp - 5._dp*z**2 - 3._dp*beta**2*(1._dp
     & - z**2)) - gvt**2*(3._dp + z**2 - beta**2*(1._dp - z**2))
     & + ((gat**2 + gvt**2)*rz*(1._dp - beta**2 - 3._dp*z**2
     & + 7._dp*beta**2*z**2 + two*beta**4*(1._dp - z**2)))/(beta**2*
     & (1._dp - beta**2)))*xI2(mt**2,mz**2,mt**2,musq,ep)
     & + two*(-4._dp*gvt**2*(1._dp + z**2) - (-3._dp*gat**2
     & + gvt**2)*f2(z,beta) + (two*(gat**2 + gvt**2)*rz*f2(z,beta))
     & /beta**2)*xI2(s,mt**2,mt**2,musq,ep) + two*s*(-(1._dp
     & + beta**2)*gvt**2*(two - beta**2*(1._dp - z**2))
     & - gat**2*(two - 3._dp*beta**2 + 5._dp*beta**2*z**2
     & + 3._dp*beta**4*(1._dp - z**2)) + (two*(gat**2
     & + gvt**2)*rz**2*f2(z,beta))/beta**2 + 4._dp*rz*(-gvt**2*
     & (1._dp + z**2) + gat**2*f2(z,beta)))*xI3(mt**2,mt**2,s,mt**2,
     & mz**2,mt**2,musq,ep)))/pi

      qa(2) =
     & (0.5_dp*alpha*gw**2*sigma0*(-two*(1._dp + z**2)
     & + (-1._dp + beta**2 + 4._dp*(-rb + rw))*s*(two - beta**2*(1._dp
     & - z**2))*db0(mt**2,mb**2,mw**2) + (8._dp*(1._dp
     & + z**2)*(-xI1(mb**2,musq,ep) + xI1(mw**2,musq,ep)))/((1._dp
     & - beta**2)*s) - (1._dp*(1._dp - 3._dp*z**2 - two*beta**4*(1._dp
     & - z**2) - 5._dp*beta**2*(1._dp + z**2) + (4._dp*(-rb + rw)*(1._dp
     & - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + two*beta**4*
     & (1._dp - z**2)))/(1._dp - beta**2))*xI2(mt**2,mb**2,mw**2,musq,
     & ep))/beta**2 + ((1._dp - 3._dp*z**2 - two*beta**4*(1._dp
     & - z**2) - 5._dp*beta**2*(1._dp + z**2) + 4._dp*(-rb + rw)*
     & f2(z,beta))*xI2(s,mb**2,mb**2,musq,ep))/beta**2 + (0.25_dp*s*
     & (1._dp - 16._dp*beta**2 + beta**4 - 3._dp*z**2 - 4._dp*beta**2*
     & z**2 - 11._dp*beta**4*z**2 - two*beta**6*(1._dp - z**2)
     & + 8._dp*rw*(1._dp - 3._dp*z**2 - two*beta**4*(1._dp - z**2)
     & - 3._dp*beta**2*(1._dp + z**2)) + 8._dp*(-(1._dp + beta**2)*rb
     & + two*(-rb + rw)**2)*f2(z,beta))*xI3(mt**2,mt**2,s,mb**2,
     & mw**2,mb**2,musq,ep))/beta**2))/pi

      qa(3) = (0.5_dp*alpha*gat**2*mt**2*sigma0*(-two*(1._dp + z**2)
     & + 4._dp*rz*s*(two - beta**2*(1._dp - z**2))*db0(mt**2,mz**2,
     & mt**2) - (8._dp*(1._dp + z**2)*(xI1(mt**2,musq,ep)
     & - xI1(mz**2,musq,ep)))/((1._dp - beta**2)*s) - (4._dp*rz*(1._dp
     & - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + two*beta**4*
     & (1._dp - z**2))*xI2(mt**2,mz**2,mt**2,musq,ep))/(beta**2*(1._dp
     & - beta**2)) + (two*(beta**2*(1._dp + z**2) + two*rz*
     & f2(z,beta))*xI2(s,mt**2,mt**2,musq,ep))/beta**2 + (4._dp*s*
     & (beta**2*(1._dp - beta**2)*rz*(1._dp - z**2)
     & + rz**2*f2(z,beta))*xI3(mt**2,mt**2,s,mt**2,mz**2,mt**2,musq,
     & ep))/beta**2))/(mz**2*pi)

      qa(4) = (0.5_dp*alpha*gw**2*sigma0*((-0.25_dp*ys*(1._dp + z**2))/rw
     & + 0.125_dp*s*((-(1._dp - beta**2)**2)/rw + 8._dp*(1._dp
     & - beta**2)*yphi - 16._dp*rb*yphi + 4._dp*ys)*(two - beta**2*
     & (1._dp - z**2))*db0(mt**2,mb**2,mw**2) + (ys*(1._dp
     & + z**2)*(-xI1(mb**2,musq,ep) + xI1(mw**2,musq,ep)))/((1._dp
     & - beta**2)*mw**2) - (0.125_dp*(-16._dp*beta**2*
     & (1._dp - beta**2)**2*yphi*(1._dp - z**2) - 16._dp*rb*yphi*(1._dp
     & - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + two*beta**4*
     & (1._dp - z**2)) + 4._dp*ys*(1._dp - beta**2 - 3._dp*z**2
     & + 7._dp*beta**2*z**2 + two*beta**4*(1._dp - z**2))
     & + ((1._dp - beta**2)**2*(1._dp - 3._dp*z**2 - two*beta**4*(1._dp
     & - z**2) + 3._dp*beta**2*(1._dp + z**2)))/rw)*xI2(mt**2,mb**2,
     & mw**2,musq,ep))/(beta**2*(1._dp - beta**2)) + (0.125_dp*(((1._dp
     & - beta**2)*(1._dp - 3._dp*z**2 - two*beta**4*(1._dp - z**2)
     & + 3._dp*beta**2*(1._dp + z**2)))/rw + 4._dp*(-two*beta**2*yphi
     & - 4._dp*rb*yphi + ys)*f2(z,beta))*xI2(s,mb**2,mb**2,musq,ep))/
     & beta**2 + (0.03125_dp*s*(8._dp*(1._dp - beta**2)**2*(1._dp
     & - 3._dp*z**2 + two*beta**2*(1._dp - z**2)) - 4._dp*(1._dp
     & - beta**2)*yphi*(1._dp - 11._dp*beta**2 - 3._dp*z**2
     & - 3._dp*beta**2*z**2 + two*beta**4*(1._dp - z**2))
     & + ((1._dp - beta**2)**2*(1._dp - 7._dp*beta**2 - 3._dp*z**2
     & + beta**2*z**2 + two*beta**4*(1._dp - z**2)))/rw
     & - 16._dp*rb*yphi*(1._dp + beta**2 - 3._dp*z**2
     & + 9._dp*beta**2*z**2 + two*beta**4*(1._dp - z**2))
     & + 16._dp*(-8._dp*rb**2 + 4._dp*rb**2*yphi + rw*ys)*f2(z,beta))*
     & xI3(mt**2,mt**2,s,mb**2,mw**2,mb**2,musq,ep))/beta**2))/pi

      qa(5) =
     & (0.5_dp*alpha*gw**2*mt**2*sigma0*(-1._dp - z**2 + two*(-1._dp
     & + beta**2 + rh)*s*(two - beta**2*(1._dp - z**2))*db0(mt**2,
     & mt**2,mh**2) + (4._dp*(1._dp + z**2)*
     & (xI1(mh**2,musq,ep) - xI1(mt**2,musq,ep)))/((1._dp - beta**2)*
     & s) - two*(two*(1._dp - beta**2)*(1._dp - z**2) + (rh*(1._dp
     & - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + two*beta**4*
     & (1._dp - z**2)))/(beta**2*(1._dp - beta**2)))*xI2(mt**2,mt**2,
     & mh**2,musq,ep) + ((beta**2*(5._dp - 3._dp*z**2 - 4._dp*beta**2*
     & (1._dp - z**2)) + two*rh*f2(z,beta))*xI2(s,mt**2,mt**2,musq,
     & ep))/beta**2 + (two*s*(3._dp*beta**2*(1._dp - beta**2)*rh*
     & (1._dp - z**2) - beta**2*(1._dp - beta**2)*(two - beta**2*
     & (1._dp - z**2)) + rh**2*f2(z,beta))*xI3(mt**2,mt**2,s,mt**2,
     & mh**2,mt**2,musq,ep))/beta**2))/(mw**2*pi)

      do ep = -1,0
      bxew(:,ep)=
     &(-0.015625_dp*as*BB(:,0)*beta*(-(one+beta*z)*xI3(zero,u,mt**2,zero,zero,mt**2,musq,ep)+(one-beta*z)*
     &xI3(zero,t,mt**2,zero,zero,mt**2,musq,ep)))/(Nc**2*pi)+(alpha*sigma0*((-two*beta**2*s*(one-z**2)*
     &(two*gaq*gat*(-rz**2+beta*z+rz*(one+beta*z))+gvq*gvt*(3._dp-beta**2+two*rz**2+two*beta*z*(one
     &+beta*z)-rz*(one+beta**2+two*beta*z)))*D6(zero,zero,mt**2,mt**2,s,u,zero,zero,mz**2,mt**2,musq,ep))/
     &((one-rz)**2*(one+beta*z))+(two*beta**2*s*(one-z**2)*(two*gaq*gat*(rz**2+beta*z-rz*(one-beta*z))
     &+gvq*gvt*(3._dp-beta**2+two*rz**2-rz*(one+beta**2-two*beta*z)-two*beta*z*(one-beta*z)))*D6(zero,zero,
     &mt**2,mt**2,s,t,zero,zero,mz**2,mt**2,musq,ep))/((one-rz)**2*(one-beta*z))-(two*(one-beta**2)*
     &(two*beta*gaq*gat+gvq*gvt*z*(one+two*beta**2-beta**2*z**2))*xI2(mt**2,zero,mt**2,musq,ep))/(beta*(one
     &-beta**2*z**2))-(two*(one-beta**2)*(two*beta*gaq*gat+gvq*gvt*z*(one+two*beta**2-beta**2*z**2))*xI2(mt**2,
     &mz**2,mt**2,musq,ep))/(beta*(one-beta**2*z**2))+(4._dp*(beta*gaq*gat+gvq*gvt*z)*xI2(s,zero,mz**2,musq,ep))
     &/beta-(two*(-gaq*gat+gvq*gvt)*(one-two*beta**2+beta**2*z**2)*xI2(u,zero,mt**2,musq,ep))/(one+beta*z)
     &+(two*(gaq*gat+gvq*gvt)*(one-two*beta**2+beta**2*z**2)*xI2(t,zero,mt**2,musq,ep))/(one-beta*z)
     &-(one*s*(two*gaq*gat*(one-beta**2+two*beta*z-beta**3*z+two*beta**2*z**2+beta**3*z**3+rz*(-3._dp
     &+3._dp*beta**2-3._dp*beta*z+beta**3*z-two*beta**2*z**2)-rz**3*(one-two*beta**2+beta**2*z**2)+rz**2*(3._dp
     &-4._dp*beta**2+beta*z-two*beta**3*z+beta**2*z**2+beta**3*z**3))+gvq*gvt*(two-beta**2+4._dp*beta*z
     &-two*beta**3*z+3._dp*beta**2*z**2-beta**4*z**2+two*beta**3*z**3+beta**4*z**4+two*rz**3*(one
     &-two*beta**2+beta**2*z**2)-rz**2*(4._dp-5._dp*beta**2-beta**4-two*beta**3*z+beta**2*z**2+beta**4*z**2
     &+two*beta**3*z**3)+beta*rz*(-4._dp*beta+beta**3-4._dp*z-two*beta**3*z**2+beta**3*z**4)))*xI3(zero,mt**2,u,zero,
     &mz**2,mt**2,musq,ep))/((one-rz)**2*(one+beta*z))+(s*(two*gaq*gat*(-one+beta**2+two*beta*z-beta**3*z
     &-two*beta**2*z**2+beta**3*z**3+rz**3*(one-two*beta**2+beta**2*z**2)+rz*(3._dp-3._dp*beta**2
     &-3._dp*beta*z+beta**3*z+two*beta**2*z**2)+rz**2*(-3._dp+4._dp*beta**2+beta*z-two*beta**3*z-beta**2*z**2
     &+beta**3*z**3))+gvq*gvt*(two-beta**2-4._dp*beta*z+two*beta**3*z+3._dp*beta**2*z**2-beta**4*z**2
     &-two*beta**3*z**3+beta**4*z**4+two*rz**3*(one-two*beta**2+beta**2*z**2)-rz**2*(4._dp-5._dp*beta**2
     &-beta**4+two*beta**3*z+beta**2*z**2+beta**4*z**2-two*beta**3*z**3)+beta*rz*(-4._dp*beta+beta**3+4._dp*z
     &-two*beta**3*z**2+beta**3*z**4)))*xI3(zero,mt**2,t,zero,mz**2,mt**2,musq,ep))/((one-rz)**2*(one-beta*z))
     &+(two*(one-beta**2)*s*(two*beta*gaq*gat*(rz**2+beta**2*z**2-rz*(one-beta**2*z**2))+gvq*gvt*z*(one
     &+two*beta**2-beta**4-beta**2*z**2+beta**4*z**2-rz**2*(-one-two*beta**2+beta**2*z**2)+rz*(-two
     &-beta**4+two*beta**2*z**2+beta**4*z**2)))*xI3(s,mt**2,mt**2,zero,mz**2,mt**2,musq,ep))/(beta*(one-rz)*(one
     &-beta**2*z**2))))/pi

      bxqcd(:,ep) =
     & (-0.015625_dp*as*BB(:,0)*beta*(-(one + beta*z)*
     & xI3(0._dp,u,mt**2,0._dp,0._dp,mt**2,musq,ep)
     & + (one - beta*z)*
     & xI3(0._dp,t,mt**2,0._dp,0._dp,mt**2,musq,ep)))
     & /(Nc**2*pi) + (alpha*sigma0*((-beta**2*s*(one - z**2)*
     & (two*beta*gaq*gat*z + gvq*gvt*(3._dp - beta**2
     & + two*beta*z*(one + beta*z)))*D6(0._dp,0._dp,mt**2,mt**2,s,
     & u,0._dp,0._dp,0._dp,mt**2,musq,ep))/(one + beta*z)
     & + (beta**2*s*(one - z**2)*(two*beta*gaq*gat*z
     & + gvq*gvt*(3._dp - beta**2 - two*beta*z*(one
     & - beta*z)))*D6(0._dp,0._dp,mt**2,mt**2,s,t,0._dp,0._dp,0._dp,
     & mt**2,musq,ep))/(one - beta*z) + (two*(one - beta**2)*
     & (-beta*gaq*gat*(one + beta**2*z**2) + gvq*gvt*z*(-one
     & - two*beta**2 + beta**2*z**2))*
     & xI2(mt**2,0._dp,mt**2,musq,ep))/(beta*(one - beta**2*z**2))
     & + (two*(beta*gaq*gat + gvq*gvt*z)*xI2(s,0._dp,0._dp,musq,ep))
     & /beta - (one*(beta*gaq*gat*(beta + z)*(one - beta*z)
     & + gvq*gvt*(one - two*beta**2 + beta**2*z**2))*xI2(u,0._dp,
     & mt**2,musq,ep))/(one + beta*z) + ((-beta*gaq*gat*(beta - z)*
     & (one + beta*z) + gvq*gvt*(one - two*beta**2
     & + beta**2*z**2))*xI2(t,0._dp,mt**2,musq,ep))/(one - beta*z)
     & + ((one - beta**2)*s*(beta*gaq*gat*(one + beta**2*z**2)
     & + gvq*gvt*z*(one + beta**2 + beta**2*(one - beta**2)*(one
     & - z**2)))*xI3(s,mt**2,mt**2,0._dp,0._dp,mt**2,musq,ep))/
     & (beta*(one - beta**2*z**2))))/(pi*(one - rz))

      end do

c check the coefficient prop to epinv
c      write(6,*) bxqcd(1:2,-1)*BB(1:2,0)/BB(1:2,-1)
c      write(6,*) as/3two/pi/Nc**2*beta/s*BB(1:2,0)*log((t-mt**2)
c     & /(u-mt**2))
c check the finite piece from epinv/epinv
c      write(6,*) bxqcd(1:2,-1)
c      write(6,*) as/3two/pi/Nc**2*beta/s*BB(1:2,-1)*log((t-mt**2)
c     & /(u-mt**2))
c      bxew(:,-1) = bxew(:,-1)*(BB(:,0)/BB(:,-1)*epinv + one)
c      bxqcd(:,-1) = bxqcd(:,-1)*(BB(:,0)/BB(:,-1)*epinv + one)

      bxew(:,-1) = bxew(:,-1)*epinv
      bxqcd(:,-1) = bxqcd(:,-1)*epinv

c--- use the value of "tevscale" (passed via anomcoup.f)
c--- as an anomalous top Yukawa coupling:
c--- g(top Y) = (tevscale) x g(SM, top Y)
c--- it therefore affects Higgs diagrams as the square
      qa(5)=qa(5)*tevscale**2

      corr = qa(1) + qa(2) + qa(3) + qa(4) + qa(5)
     & + bxew(:,0) + bxew(:,-1) + bxqcd(:,0) + bxqcd(:,-1)
      corr = corr + ini_corr

      corr = corr/born

      end subroutine qqbQQb_ew_oneloop


      function f1(x)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'cplx.h'
      real(dp):: x,cli2,f1

      f1 = real(one + two*(
     & + (one+log(cplx1(x)))*(two*x+3._dp)
     & - two*(one+x)**2*(cli2(cplx1(one+one/x))
     & -pi**2/6._dp)),dp)

      end function f1


      function f2(z,beta)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: z,beta,f2

      f2 = one - 3._dp*z**2 - two*beta**2*(one - z**2)

      end function f2


      function D6(p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,
     & m2sq,m3sq,m4sq,musq,ep)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,
     & m4sq,musq,xI3,xI4,p1dp2,p1dp3,p2dp3,m1,m2,m3,m4,D6
      integer ep
      external xI3,xI4

      p1dp2 = (s12 - p1sq - p2sq)/two
      p2dp3 = (s23 - p2sq - p3sq)/two
      p1dp3 = (p2sq + p4sq - s12 -s23)/two
      m1 = sqrt(m1sq)
      m2 = sqrt(m2sq)
      m3 = sqrt(m3sq)
      m4 = sqrt(m4sq)

c23456789012345678901234567890123456789012345678901234567890123456789012
      D6=(0.25_dp*(-((-m1sq+m2sq)*p1dp3*p2sq+p1dp2*((m1sq-m2sq-p1sq)
     &*p2dp3+p1dp3*(m2sq-m3sq+p2sq))-p1dp2**2*(m3sq-m4sq+two*p2dp3+p3sq)
     &+p1sq*((-m2sq+m3sq)*p2dp3+p2sq*(m3sq-m4sq+p1dp3+p2dp3+p3sq)))
     &*xI3(p1sq,p2sq,s12,m1sq,m2sq,m3sq,musq,ep)+(p1dp3**2*(m2sq-m3sq
     & +p2sq)-p1dp3*(m1sq-m2sq-p1sq)*(p2dp3+p2sq)-p1dp2**2*(m3sq-m4sq
     & +two*p2dp3+p3sq)+p1sq*(two*p2dp3**2+(m3sq-m4sq+p2dp3)*p2sq-p2dp3
     & *(m2sq-two*m3sq+m4sq-p3sq)+(-m2sq+m3sq)*p3sq)+p1dp2*(p1dp3*(m2sq
     & -two*m3sq+m4sq-two*p2dp3+p2sq-p3sq)+(m1sq-m2sq-p1sq)*(p2dp3+p3sq)
     & ))*xI3(p1sq,s23,p4sq,m1sq,m2sq,m4sq,musq,ep)-((m1sq-m2sq-p1sq)
     & *(p2dp3**2-p2sq*p3sq)+p1dp2*(-two*p2dp3**2+(m2sq-m3sq+p2sq)*p3sq
     & -p2dp3*(m3sq-m4sq+p3sq))+p1dp3*((-m2sq+m3sq)*p2dp3+p2sq*(m3sq
     & -m4sq+p2dp3+p3sq)))*xI3(p2sq,p3sq,s23,m2sq,m3sq,m4sq,musq,ep)
     & +(-m3sq*p1sq*p2dp3+m4sq*p1sq*p2dp3+m1sq*p2dp3**2-m2sq*p2dp3**2
     & -p1sq*p2dp3**2+p1dp3**2*(-m2sq+m3sq+p2sq)+two*p1dp2**2*p3sq+m2sq
     & *p1sq*p3sq-m3sq*p1sq*p3sq-p1sq*p2dp3*p3sq-m1sq*p2sq*p3sq+m2sq
     & *p2sq*p3sq+p1dp2*(-two*p2dp3**2+(-m1sq+two*m2sq-m3sq+p1sq+p2sq)
     & *p3sq-p2dp3*(m3sq-m4sq+p3sq)+p1dp3*(m3sq-m4sq-two*p2dp3+p3sq))
     & +p1dp3*((m1sq-two*m2sq+m3sq-p1sq)*p2dp3+p2sq*(m3sq-m4sq+p2dp3
     & +p3sq)))*xI3(s12,p3sq,p4sq,m1sq,m3sq,m4sq,musq,ep)+(-two*m2sq
     & *m3sq*p1sq*p2dp3+two*m3**4*p1sq*p2dp3+two*m2sq*m4sq*p1sq*p2dp3
     & -two*m3sq*m4sq*p1sq*p2dp3-m1**4*p2dp3**2+two*m1sq*m2sq*p2dp3**2
     & -m2**4*p2dp3**2+two*m1sq*p1sq*p2dp3**2-two*m2sq*p1sq*p2dp3**2
     & +4._dp*m3sq*p1sq*p2dp3**2-p1sq**2*p2dp3**2+m3**4*p1sq*p2sq-two
     & *m3sq*m4sq*p1sq*p2sq+m4**4*p1sq*p2sq+two*m3sq*p1sq*p2dp3*p2sq-two
     & *m4sq*p1sq*p2dp3*p2sq-p1dp3**2*((m2sq-m3sq)**2-two*(m2sq+m3sq)
     & *p2sq+p2sq**2)+m2**4*p1sq*p3sq-two*m2sq*m3sq*p1sq*p3sq+m3**4*p1sq
     & *p3sq-two*m2sq*p1sq*p2dp3*p3sq+two*m3sq*p1sq*p2dp3*p3sq+m1**4
     & *p2sq*p3sq-two*m1sq*m2sq*p2sq*p3sq+m2**4*p2sq*p3sq-two*m1sq*p1sq
     & *p2sq*p3sq-two*m4sq*p1sq*p2sq*p3sq+p1sq**2*p2sq*p3sq+two*p1sq
     & *p2dp3*p2sq*p3sq+p1sq*p2sq**2*p3sq+p1sq*p2sq*p3sq**2-p1dp2**2
     & *((m3sq-m4sq)**2+4._dp*p2dp3**2-two*(two*m2sq-m3sq+m4sq)*p3sq
     & +p3sq**2+4._dp*p2dp3*(m3sq-m4sq+p3sq))-two*p1dp3*(m1sq-m2sq-p1sq)
     & *((-m2sq+m3sq)*p2dp3+p2sq*(m3sq-m4sq+p2dp3+p3sq))+two*p1dp2
     & *((m1sq-m2sq-p1sq)*(two*p2dp3**2-(m2sq-m3sq+p2sq)*p3sq+p2dp3
     & *(m3sq-m4sq+p3sq))+p1dp3*(-two*(m2sq+m3sq)*p2dp3+(m2sq-m3sq)
     & *(m3sq-m4sq+p3sq)+p2sq*(m3sq-m4sq+two*p2dp3+p3sq))))
     & *xI4(p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq,musq,ep)))
     & /(-two*p1dp2*p1dp3*p2dp3+p1dp3**2*p2sq+p1dp2**2*p3sq
     &+p1sq*(p2dp3**2-p2sq*p3sq))

      D6 = -two*pi*D6

      end function D6


