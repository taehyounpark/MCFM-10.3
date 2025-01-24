!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_higgs_scalar(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     f(-p1) + f(-p2) --> H + f(p5)
c                   |
c                   --> b(p3)+bbar(p4)

c--all momenta incoming

c--- Matrix elements are taken from:
c--- R.~K.~Ellis, I.~Hinchliffe, M.~Soldate and J.~J.~van der Bij,
c--- %``Higgs Decay To Tau+ Tau-: A Possible Signature Of Intermediate
c--- % Mass Higgs Bosons At The SSC,''
c--- Nucl.\ Phys.\ B {\bf 297}, 221 (1988).
      include 'constants.f'
      include 'masses.f'
      include 'nf.f'
      include 'hdecaymode.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'scalarselect.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,qq,hdecay
      real(dp):: origmbsq,dotvec,s34,qqTotal,ggTotal
      real(dp):: msqgamgam

      s34=dotvec(p(3,:)+p(4,:),p(3,:)+p(4,:))

c************ stop variables
c      mass_gap=800
c      mscalar1=mt

c      theta = pi/4d0
c      tanBeta = 10
c************

c   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
        call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
        call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
        hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in qqb_higgs'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      origmbsq=mbsq
      mbsq=mt**2

      gg=+avegg*ggTotal(s(1,2),s(1,5),s(2,5))*hdecay
      qq=+aveqq*qqTotal(s(1,2),s(1,5),s(2,5))*hdecay
      qg=-aveqg*qqTotal(s(1,5),s(1,2),s(2,5))*hdecay
      gq=-aveqg*qqTotal(s(2,5),s(1,5),s(1,2))*hdecay

      mbsq=origmbsq

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp

      if ((j== 0) .or. (k==0)) then
         if ((j== 0) .and. (k==0)) then
            msq(j,k)=gg
         elseif ((j==0).and.(k /= 0)) then
            msq(j,k)=gq
         elseif ((j /= 0).and.(k==0)) then
            msq(j,k)=qg
         endif
      elseif ((j==-k).and. (j /= 0)) then
         msq(j,k)=qq
      endif

      enddo
      enddo

      return
      end



c ********************************************* !
c ****************** Totals ******************* !
c ********************************************* !

      function ggTotal(s,t,u)
      implicit none
      include 'types.f'
      real(dp):: ggTotal
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      real(dp):: s,t,u
      real(dp):: ppp,ppm

      ggTotal = 2d0*ppp(s,t,u)
     &           + 2d0*ppm(s,t,u)
     &           + 2d0*ppp(t,s,u)
     &           + 2d0*ppp(u,t,s)

      ggTotal = 3d0*as**3*gwsq/32d0/pi/wmass**2*ggTotal

      return
      end

      function qqTotal(s,t,u)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      real(dp):: qqTotal
      real(dp):: s,t,u,s1
      complex(dp):: qqscalar,bgqqg

      s1=s-hmass**2

      qqTotal = abs(qqscalar(1,s,t,u)
     &           + qqscalar(2,s,t,u)
     &           + bgqqg(s,t,u))**2

      qqTotal = (t**2+u**2)/s/s1**2*qqTotal

      qqTotal = 4*as**3*gwsq/pi/wmass**2*qqTotal

      return
      end

      function ppp(s,t,u)
      implicit none
      include 'types.f'
      real(dp):: ppp
      complex(dp)::bgppp,pppScalar
      real(dp)::s,t,u

      ppp = abs(bgppp(s,t,u)+pppScalar(1,s,t,u)+pppScalar(2,s,t,u))**2


      return
      end

      function ppm(s,t,u)
      implicit none
      include 'types.f'
      real(dp):: ppm
      complex(dp)::bgppm,ppmScalar
      real(dp)::s,t,u

      ppm = abs(bgppm(s,t,u)+ppmScalar(1,s,t,u)+ppmScalar(2,s,t,u))**2

      return
      end



c ********************************************* !
c **************** Fermionic ****************** !
c ********************************************* !

      function bgppp(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: bgppp
      include 'masses.f'
c---Matrix element squared eq. A.15
      real(dp):: delta,s1,t1,u1,s,t,u,mtsq,hmass2
      complex(dp):: BoxInt,C1,B1,E

      delta = (s*t*u/8)**0.5
      hmass2 = hmass**2
      s1 = s-hmass2
      t1 = t-hmass2
      u1 = u-hmass2

      mtsq=mt**2

      bgppp= -64*(1/u/t+1/t/t1+1/u/u1)
     &        -64/s*((2*s+t)/u1**2*B1(u,mtsq)+(2*s+u)/t1**2*B1(t,mtsq))
     &        -16*(s-4*mtsq)/s/t/u * ( s1*C1(s,mtsq) + (u-s)*C1(t,mtsq) + (t-s)*C1(u,mtsq) )
     &        -128*mtsq*(1/t/t1*C1(t,mtsq) + 1/u/u1*C1(u,mtsq))
     &        +64*mtsq/s*BoxInt(u,t,mtsq)
     &        +8*(s-4*mtsq)/s/t/u*(s*t*BoxInt(s,t,mtsq)+ u*s*BoxInt(u,s,mtsq)-u*t*BoxInt(u,t,mtsq))
     &        -32/s**2*E(u,t,mtsq)

      bgppp = mtsq*delta*bgppp

      return
      end

      function bgppm(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: bgppm
      include 'masses.f'
c---Matrix element squared eq. A.16
      real(dp):: s,t,u,delta,s1,t1,u1,mtsq,hmass2
      complex(dp):: BoxInt,C1

      delta = (s*t*u/8)**0.5
      mtsq=mt**2
      hmass2 = hmass**2

      s1 = s-hmass2
      t1 = t-hmass2
      u1 = u-hmass2

      bgppm= 64*hmass2/s/t/u
     &        +16*(hmass2-4*mtsq)/s/t/u*( s1*C1(s,mtsq) + u1*C1(u,mtsq) + t1*C1(t,mtsq) )
     &        -8*(hmass2-4*mtsq)/s/t/u*( s*t*BoxInt(s,t,mtsq)+u*s*BoxInt(u,s,mtsq)+u*t*BoxInt(u,t,mtsq) )

      bgppm = mtsq*delta*bgppm
      return
      end

      function bgqqg(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: bgqqg

      include 'masses.f'
c---Matrix element squared eq. A.15
      real(dp):: s1,s,t,u,hmass2
      complex(dp):: C1,B1

      hmass2 = hmass**2

      s1 = s-hmass2

      bgqqg = mt**2*(2d0+2d0*s/s1*B1(s,mt**2)+(4d0*mt**2-u-t)*C1(s,mt**2))

      return
      end















c ********************************************* !
c ****************** Scalar ******************* !
c ********************************************* !

      function pppscalar(i,s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: pppscalar
      include 'masses.f'
      include 'ewcouple.f'
      include 'squark.f'
      integer i
      real(dp)::s,t,u,s1,t1,u1,delta
      real(dp)::msq
      complex(dp)::B1,C1,BoxInt,E

      delta = (s*t*u/8)**0.5
      s1 = s-hmass**2
      t1 = t-hmass**2
      u1 = u-hmass**2

      msq=mloopsq(i+1)

      pppscalar = 16*(1/t/u+1/t1/t+1/u1/u)
     &             + 16/s*( B1(t,msq)*(2*s+u)/t1**2 + B1(u,msq)*(2*s+t)/u1**2 )
     &             + 32*msq*( C1(t,msq)/t/t1 + C1(u,msq)/u/u1 )
     &             - 16*msq/s/t/u*( s1*C1(s,msq)+(u-s)*C1(t,msq)+(t-s)*C1(u,msq) )
     &             + 8*msq/s/t/u*( s*t*BoxInt(s,t,msq)+s*u*BoxInt(s,u,msq)-t*u*BoxInt(t,u,msq) )
     &             - 16*msq/s*BoxInt(t,u,msq) +8/s**2*E(t,u,msq)

c      pppscalar = g(i,theta,tanBeta,deltaMsq)*vevsq**0.5*delta*pppscalar
      pppscalar = -4d0*gloop(i+1)*vevsq**0.5*delta*pppscalar

      return
      end

      function ppmscalar(i,s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ppmscalar
      include 'masses.f'
      include 'ewcouple.f'
      include 'squark.f'
      integer i
      real(dp)::s,t,u,s1,t1,u1,delta
      real(dp)::msq
      complex(dp)::C1,BoxInt

      delta = (s*t*u/8)**0.5
      s1 = s-hmass**2
      t1 = t-hmass**2
      u1 = u-hmass**2

      msq=mloopsq(i+1)

      ppmscalar = - 16*hmass**2/s/t/u
     &            + 16*msq/s/t/u*( s1*C1(s,msq) + t1*C1(t,msq) + u1*C1(u,msq) )
     &            - 8*msq/s/t/u*( s*t*BoxInt(s,t,msq) +s*u*BoxInt(s,u,msq) + t*u*BoxInt(t,u,msq) )

c      ppmscalar = g(i,theta,tanBeta,deltaMsq)*vevsq**0.5*delta*ppmscalar
      ppmscalar = -4d0*gloop(i+1)*vevsq**0.5*delta*ppmscalar

      return
      end

      function qqscalar(i,s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: qqscalar
      include 'masses.f'
      include 'ewcouple.f'
      include 'squark.f'
      real(dp):: s,t,u,s1,msq

      integer i
      complex(dp):: C1,B1

      s1 = s-hmass**2

      msq=mloopsq(i+1)

c      qqscalar = -g(i,theta,tanBeta,deltaMsq)*vevsq**0.5
c     &            * (1/2d0+msq*C1(s,msq)+s/2/s1*B1(s,msq))
      qqscalar = 4d0*gloop(i+1)*vevsq**0.5
     &          * (1/2d0+msq*C1(s,msq)+s/2/s1*B1(s,msq))


      return
      end

      function g(i,theta,tanBeta,deltaMsq)
      implicit none
      include 'types.f'
      real(dp):: g
      include 'masses.f'
      include 'ewcouple.f'
      integer i
      real(dp)::alpha1,alpha2,theta,tanBeta,deltaMsq
      real(dp)::cTheta,sTheta,c2Beta


      c2Beta = cos(2*atan(tanBeta))

      cTheta = cos(theta)
      sTheta = sin(theta)

      alpha1 = zmass**2/mt**2*c2Beta*(1d0-4d0/3d0*xw)
      alpha2 = 4d0/3d0*zmass**2/mt**2*c2Beta*xw

      if (i==1) then
      g = mt**2/sqrt(vevsq)
     & *(alpha1*cTheta**2+alpha2*sTheta**2
     & +2._dp-deltaMsq/2/mt**2*sin(2*theta)**2)
      elseif (i==2) then
      g = mt**2/sqrt(vevsq)
     & *(alpha1*sTheta**2+alpha2*cTheta**2
     & +2._dp+deltaMsq/2/mt**2*sin(2*theta)**2)
      else
      g=0._dp
      stop 'function g=0'
      endif
      return
      end function g












c ********************************************* !
c **************** Integrals ****************** !
c ********************************************* !


      function E(u,t,masssq)
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      ! use iso_c_binding
      ! use iso_fortran_env
      implicit none
      include 'types.f'
      complex(dp):: E,C1,BoxInt
      include 'masses.f'
      include 'scale.f'
c--- eq. B.12
      real(dp):: t,u,t1,u1,hmass2,masssq

      hmass2 = hmass**2

      t1 = t-hmass2
      u1 = u-hmass2

      E =u*loopI3(0d0,0d0,u,masssq,masssq,masssq,musq,0)
     &     + t*loopI3(0d0,0d0,t,masssq,masssq,masssq,musq,0)
     &     + u1*C1(u,masssq) + t1*C1(t,masssq)
     &     - u*t*BoxInt(u,t,masssq)

      return
      end

      function C1(s,masssq)
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      ! use iso_c_binding
      ! use iso_fortran_env
      implicit none
      include 'types.f'
      complex(dp):: C1
      include 'masses.f'
      include 'scale.f'
c--- eq. B.7
      real(dp):: s,s1,masssq,hmass2
      hmass2 = hmass**2

      s1=s-hmass2

      C1 = (s*loopI3(0d0,0d0,s,masssq,masssq,masssq,musq,0) - hmass2*loopI3(0d0,0d0,hmass2,masssq,masssq,masssq,musq,0))/s1
      return
      end

      function B1(s,masssq)
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      ! use iso_c_binding
      ! use iso_fortran_env
      implicit none
      include 'types.f'
      complex(dp):: B1
      include 'masses.f'
      include 'scale.f'
c--- eq. B.5
      real(dp):: s,hmass2,masssq

      hmass2 = hmass**2


      B1 = loopI2(s,masssq,masssq,musq,0) -loopI2(hmass2,masssq,masssq,musq,0)

      return
      end

      function BoxInt(s,t,masssq)
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      ! use iso_c_binding
      ! use iso_fortran_env
      implicit none
      include 'types.f'
      complex(dp):: BoxInt
      include 'masses.f'
      include 'scale.f'
      real(dp):: s,t,masssq,hmass2

      hmass2 = hmass**2

      BoxInt = loopI4(0d0,0d0,0d0,hmass2,s,t,masssq,masssq,masssq,masssq,musq,0)

      return
      end
