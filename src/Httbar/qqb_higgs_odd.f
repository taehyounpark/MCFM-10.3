!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_higgs_odd(p,msq)
      implicit none
      include 'types.f'

c---Matrix element squared averaged over initial colors and spins
c     f(-p1) + f(-p2) --> A + f(p5)
c                         |
c                         --> b(p3)+bbar(p4)

c   Here A is pseudoscalar (CP-odd) Higgs
c   Note that the approximate forms (for large mbsq) are exactly
c   the same as in the scalar case, apart from an overall factor
c--all momenta incoming
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'hdecaymode.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,qq,hdecay
      real(dp):: ehsvm3_odd,ehsvm4_odd,origmbsq,s34,msqgamgam
      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

      msq(:,:)=zip

c   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
          write(6,*) 'Unimplemented process in qqb_higgs_odd'
          stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      origmbsq=mbsq
      mbsq=mt**2
      gg=+avegg*ehsvm3_odd(s(1,2),s(1,5),s(2,5))*hdecay
      qq=+aveqq*ehsvm4_odd(s(1,2),s(1,5),s(2,5))*hdecay
      qg=-aveqg*ehsvm4_odd(s(1,5),s(1,2),s(2,5))*hdecay
      gq=-aveqg*ehsvm4_odd(s(2,5),s(1,5),s(1,2))*hdecay
      mbsq=origmbsq

      do j=-nf,nf
      do k=-nf,nf

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


      function ehsvm3_odd(s,t,u)
      implicit none
      include 'types.f'
      real(dp):: ehsvm3_odd

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
c---Matrix element squared Eqn 2.2 of EHSV
      complex(dp):: ehsva2_odd,ehsva4_odd
      real(dp):: hmass2,s,t,u
      logical:: approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq

      hmass2=s+t+u
      if (approx) then
      ehsvm3_odd=gwsq/pi*as**3*xn*V/9._dp*(
     &        hmass2**4+s**4+t**4+u**4)/s/t/u/wmass**2
      else
      ehsvm3_odd=
     & abs(ehsva2_odd(s,t,u))**2+abs(ehsva2_odd(u,s,t))**2
     &+abs(ehsva2_odd(t,u,s))**2+abs(ehsva4_odd(s,t,u))**2
      ehsvm3_odd=gwsq/pi*as**3*xn*V*hmass2**4/(s*t*u*wmass**2)*ehsvm3_odd
      endif
c--- effective Lagrangian increased by factor of (3/2) wrt. scalar case,
c--- hence that factor squared in the rate
      ehsvm3_odd=ehsvm3_odd*(3._dp/2._dp)**2

      return
      end


      function ehsvm4_odd(s,t,u)
      implicit none
      include 'types.f'
      real(dp):: ehsvm4_odd

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
c---Matrix element squared Eqn 2.6 of EHSV
      complex(dp):: ehsva5_odd
      real(dp):: hmass2,s,t,u
      logical:: approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq

      hmass2=s+t+u
      if (approx) then
      ehsvm4_odd=gwsq/pi*as**3*V/18._dp*(u**2+t**2)/s/wmass**2
      else
      ehsvm4_odd=abs(ehsva5_odd(s,t,u))**2
      ehsvm4_odd=gwsq/(4._dp*pi)*as**3*V/2._dp*(u**2+t**2)/(s*wmass**2)
     & *hmass2**2/(u+t)**2*ehsvm4_odd
      endif
c--- effective Lagrangian increased by factor of (3/2) wrt. scalar case,
c--- hence that factor squared in the rate
      ehsvm4_odd=ehsvm4_odd*(3._dp/2._dp)**2

      return
      end



