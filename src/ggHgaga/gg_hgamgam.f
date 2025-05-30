!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_hgamgam(p,msq)
      implicit none
      include 'types.f'

c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
c----averaged over initial colours and spins
c    g(-p1)+g(-p2)-->H --> gamma(p3) + gamma(p4)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),s12
      real(dp):: decay,gg,Asq,msqgamgam
c---set msq=0 to initialize
      msq(:,:)=0._dp

      s12=2._dp*(p(1,4)*p(2,4)-p(1,1)*p(2,1)-p(1,2)*p(2,2)-p(1,3)*p(2,3))
      decay=msqgamgam(hmass)/((s12-hmass**2)**2+(hmass*hwidth)**2)
      Asq=(as/(3._dp*pi))**2/vevsq
      gg=0.5_dp*V*Asq*s12**2

c---calculate propagators
      msq(0,0)=avegg*gg*decay

      return
      end


