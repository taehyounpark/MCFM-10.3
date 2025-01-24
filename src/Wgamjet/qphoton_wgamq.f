!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qphoton_wgamq(p,msq)
      implicit none
      include 'types.f'

c-----Author Keith Ellis, May 2020
c----Matrix element for W gam production
c----averaged over initial colours and spins
c For nwz=+1
c     u(-p1)+photon(-p2)-->W^+(n(p3)+e^+(p4)) + gamma(p5) + d(p6)
c For nwz=-1
c     ubar(-p1)+photon(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5) + dbar(p6)
c---
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ckm.f'
      include 'zprods_com.f'
      include 'nwz.f'
      include 'zcouple_cms.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: qgsq,gqsq,agsq,gasq
      complex(dp)::gq(2,2),qg(2,2),ga(2,2),ag(2,2)

      call spinoru(6,p,za,zb)
      fac=xn*abs((zesq/zxw)**2*zesq**2)

      if (nwz == -1) then
        call a6Wpgamgam(6,1,4,3,5,2,zb,za,qg)
        call a6Wpgamgam(2,6,4,3,5,1,zb,za,ga)
        call a6Wpgamgam(6,2,4,3,5,1,zb,za,gq)
        call a6Wpgamgam(1,6,4,3,5,2,zb,za,ag)
c These are included as part of the Wga_ew process
c        call a6Wpgamgam(2,1,4,3,5,6,zb,za,qa)
c        call a6Wpgamgam(1,2,4,3,5,6,zb,za,aq)
      elseif (nwz == +1) then
        call a6Wpgamgam(1,6,3,4,5,2,za,zb,qg)
        call a6Wpgamgam(6,2,3,4,5,1,za,zb,ga)
        call a6Wpgamgam(2,6,3,4,5,1,za,zb,gq)
        call a6Wpgamgam(6,1,3,4,5,2,za,zb,ag)
c These are included as part of the Wga_ew process
c        call a6Wpgamgam(1,2,3,4,5,6,za,zb,qa)
c        call a6Wpgamgam(2,1,3,4,5,6,za,zb,aq)
      endif

      qgsq=spinave/xn*fac*(cdabs(qg(1,1))**2+cdabs(qg(1,2))**2
     &                    +cdabs(qg(2,1))**2+cdabs(qg(2,2))**2)
      gasq=spinave/xn*fac*(cdabs(ga(1,1))**2+cdabs(ga(1,2))**2
     &                    +cdabs(ga(2,1))**2+cdabs(ga(2,2))**2)
      gqsq=spinave/xn*fac*(cdabs(gq(1,1))**2+cdabs(gq(1,2))**2
     &                    +cdabs(gq(2,1))**2+cdabs(gq(2,2))**2)
      agsq=spinave/xn*fac*(cdabs(ag(1,1))**2+cdabs(ag(1,2))**2
     &                    +cdabs(ag(2,1))**2+cdabs(ag(2,2))**2)
c These are included as part of the Wga_ew process
c      qasq=half*aveqq*fac*(cdabs(qa(1,1))**2+cdabs(qa(1,2))**2
c     &                    +cdabs(qa(2,1))**2+cdabs(qa(2,2))**2)
c      aqsq=half*aveqq*fac*(cdabs(aq(1,1))**2+cdabs(aq(1,2))**2
c     &                    +cdabs(aq(2,1))**2+cdabs(aq(2,2))**2)

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize and abuse 0 to mean photon instead of gluon
      msq(j,k)=0._dp
c          if ((j > 0) .and. (k < 0)) then
c            msq(j,k)=Vsq(j,k)*qasq*two    ! note factor of two
c          elseif ((j < 0) .and. (k > 0)) then
c            msq(j,k)=Vsq(j,k)*aqsq*two    ! note factor of two
          if     ((j == 0) .and. (k > 0)) then
            msq(j,k)=Vsum(k)*gqsq
          elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=Vsum(k)*gasq
          elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*qgsq
         elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*agsq
          endif
      enddo
      enddo

      return
      end


