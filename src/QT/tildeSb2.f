!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tildeSb2(i,tSb2)
c     first argument is gg or qq
c     second argument is power of Lb
c     1909.00811v2, Eq 3.16
      implicit none
      include 'types.f'
      include 'qtconstants.f'
      include 'first.f'
      include 'Lnu.f'
      include 'mpicommon.f'
      real(dp)::tSb2(0:4)
      integer i
      if (first) then
!$omp master
        if (rank == 0) then
      write(6,*) 'tildeSb2:Lnu',Lnu
      first=.false.
        endif
!$omp end master
      endif
      tSb2(4)=Gamma0(i)**2/8._dp
      tSb2(3)=-Gamma0(i)*(Lnu*Gamma0(i)+beta0/3._dp
     & +0.25_dp*(tgammaS0(i)+tgamman0(i)))
      tSb2(2)=Lnu**2*2*Gamma0(i)**2
     & +Lnu*Gamma0(i)*(beta0+tgammaS0(i)+1.5_dp*tgamman0(i))
     & +beta0*(tgammaS0(i)/4._dp+tgamman0(i)/2._dp)
     & +1/8._dp*(tgammaS0(i)+tgamman0(i))**2
     & -Gamma1(i)/2._dp-Gamma0(i)/2._dp*tildes1(i)
      tSb2(1)=-Lnu**2*2*Gamma0(i)*tgamman0(i)
     & +Lnu*(-(beta0+0.5d0*(tgammaS0(i)+tgamman0(i)))*tgamman0(i)
     & +2*Gamma1(i)+2*Gamma0(i)*tildes1(i))
     & +0.5d0*(tgammaS1(i)+tgamman1(i))
     & +(beta0+0.5d0*(tgammaS0(i)+tgamman0(i)))*tildes1(i)
      tSb2(0)=Lnu**2*tgamman0(i)**2/2._dp
     & -Lnu*(tgamman1(i)+tgamman0(i)*tildes1(i))+tildes2(i)
      return
      end
