!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_zh_gaga(p,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     December, 1998.                                                  *
c  Matrix element squared averaged over initial colors and spins       *
c     q(-p1)+qbar(-p2) -->  H  + Z                                     *
c                           |    |                                     *
c                           |     ->fermion(p3)+antifermion(p4)        *
c                           |                                          *
c                            ---> gamma(p5)+gamma(p6)                  *
c***********************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      include 'zcouple_cms.f'
      include 'ewinput.f'
      integer:: j,k
      real(dp):: p(mxpart,4)
      real(dp):: s,prop,fac,q1423,q2413,s56
      real(dp):: msq(-nf:nf,-nf:nf),hdecay
      real(dp):: msqhgamgam
      complex(dp):: v2(2)

      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      v2(1)=zl1
      v2(2)=zr1

      msq(:,:)=0._dp
      if (hdecaymode=='none') then
         s56=s(5,5)
      else
         s56=s(5,6)
      endif

c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)

      if (ewscheme < 4) then
        fac=xn*4._dp*(xw/(one-xw))**2*gwsq**3*wmass**2/prop
      else
c-- note: in CMS, wmass**2 -> sqrt(wmass**4+(wmass*wwidth)**2)
        fac=xn*4._dp*abs(zxw/(cone-zxw))**2*gwsq**3*sqrt(wmass**4+(wmass*wwidth)**2)/prop
      endif

      if (hdecaymode == 'none') then
         hdecay = one
      else
         hdecay=msqhgamgam(s56)
         hdecay=hdecay/((s56-hmass**2)**2+(hmass*hwidth)**2)
      endif
      fac=fac*hdecay

c-- Old form of this matrix element (modified to facilitate extension
c--- to H->WW decay)
c      spinave only (color factors cancel)
c       fac=two*spinave*gw**8*(xw/(1._dp-xw))**2*mbsq*(s56-4._dp*mb**2)/prop


      q1423=aveqq*fac*s(1,4)*s(2,3)
      q2413=aveqq*fac*s(2,4)*s(1,3)

      do j=-nf,nf
      if (j == 0) go to 40
      k=-j
      if ((j > 0) .and. (k < 0)) then
      msq(j,k)=
     &  +(abs(zl(j)*v2(1))**2+abs(zr(j)*v2(2))**2)*q1423
     &  +(abs(zl(j)*v2(2))**2+abs(zr(j)*v2(1))**2)*q2413
      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=
     &  +(abs(zl(k)*v2(1))**2+abs(zr(k)*v2(2))**2)*q2413
     &  +(abs(zl(k)*v2(2))**2+abs(zr(k)*v2(1))**2)*q1423
      endif
 40   continue
      enddo
      return
      end

c + L1(j)*L1(j)*L2(j)*L2(j)*gw^8*[1/4/XN]*mb^2*[s56-4*mb^2]
c  * ( 2*zprop1^-1*zprop2^-1*hprop^-1*sinw^4*s14*s23*cos^-4 )

c + R1(j)*R1(j)*R2(j)*R2(j)*gw^8*[1/4/XN]*mb^2*[s56-4*mb^2]
c  * ( 2*zprop1^-1*zprop2^-1*hprop^-1*sinw^4*s14*s23*cos^-4 )

c + L1(j)*L1(j)*R2(j)*R2(j)*gw^8*[1/4/XN]*mb^2*[s56-4*mb^2]
c  * ( 2*zprop1^-1*zprop2^-1*hprop^-1*sinw^4*s13*s24*cos^-4 )

c + R1(j)*R1(j)*L2(j)*L2(j)*gw^8*[1/4/XN]*mb^2*[s56-4*mb^2]
c  * ( 2*zprop1^-1*zprop2^-1*hprop^-1*sinw^4*s13*s24*cos^-4 )


