!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wgamg_gam(p,msq)
      implicit none
      include 'types.f'

c-----Author Keith Ellis, May 2020
c----Matrix element for W gam production
c----averaged over initial colours and spins
c For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4)) + gamma(p5) + g(p6) + gamma(p7)
c For nwz=-1
c     ubar(-p1)+d(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5) + g(p6) + gamma(p7)
c---
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'zprods_com.f'
      include 'nwz.f'
      include 'zcouple_cms.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: aqsq,qasq,qgsq,gqsq,agsq,gasq
      complex(dp)::qa(2,2,2),aq(2,2,2),gq(2,2,2),qg(2,2,2),ga(2,2,2),ag(2,2,2)

      call spinoru(7,p,za,zb)
      fac=V/2._dp*gsq*abs((zesq/zxw)**2*zesq**2)

      if (nwz == -1) then
        call a7Wpgamgamg(2,1,4,3,5,7,6,zb,za,qa)
        call a7Wpgamgamg(1,2,4,3,5,7,6,zb,za,aq)
        call a7Wpgamgamg(6,1,4,3,5,7,2,zb,za,qg)
        call a7Wpgamgamg(2,6,4,3,5,7,1,zb,za,ga)
        call a7Wpgamgamg(6,2,4,3,5,7,1,zb,za,gq)
        call a7Wpgamgamg(1,6,4,3,5,7,2,zb,za,ag)
      elseif (nwz == +1) then
        call a7Wpgamgamg(1,2,3,4,5,7,6,za,zb,qa)
        call a7Wpgamgamg(2,1,3,4,5,7,6,za,zb,aq)
        call a7Wpgamgamg(1,6,3,4,5,7,2,za,zb,qg)
        call a7Wpgamgamg(6,2,3,4,5,7,1,za,zb,ga)
        call a7Wpgamgamg(2,6,3,4,5,7,1,za,zb,gq)
        call a7Wpgamgamg(6,1,3,4,5,7,2,za,zb,ag)
      endif

      qasq=aveqq*fac*(cdabs(qa(1,1,1))**2+cdabs(qa(1,1,2))**2
     &               +cdabs(qa(1,2,1))**2+cdabs(qa(1,2,2))**2
     &               +cdabs(qa(2,1,1))**2+cdabs(qa(2,1,2))**2
     &               +cdabs(qa(2,2,1))**2+cdabs(qa(2,2,2))**2)
      aqsq=aveqq*fac*(cdabs(aq(1,1,1))**2+cdabs(aq(1,1,2))**2
     &               +cdabs(aq(1,2,1))**2+cdabs(aq(1,2,2))**2
     &               +cdabs(aq(2,1,1))**2+cdabs(aq(2,1,2))**2
     &               +cdabs(aq(2,2,1))**2+cdabs(aq(2,2,2))**2)
      qgsq=aveqg*fac*(cdabs(qg(1,1,1))**2+cdabs(qg(1,1,2))**2
     &               +cdabs(qg(1,2,1))**2+cdabs(qg(1,2,2))**2
     &               +cdabs(qg(2,1,1))**2+cdabs(qg(2,1,2))**2
     &               +cdabs(qg(2,2,1))**2+cdabs(qg(2,2,2))**2)
      gasq=aveqg*fac*(cdabs(ga(1,1,1))**2+cdabs(ga(1,1,2))**2
     &               +cdabs(ga(1,2,1))**2+cdabs(ga(1,2,2))**2
     &               +cdabs(ga(2,1,1))**2+cdabs(ga(2,1,2))**2
     &               +cdabs(ga(2,2,1))**2+cdabs(ga(2,2,2))**2)
      gqsq=aveqg*fac*(cdabs(gq(1,1,1))**2+cdabs(gq(1,1,2))**2
     &               +cdabs(gq(1,2,1))**2+cdabs(gq(1,2,2))**2
     &               +cdabs(gq(2,1,1))**2+cdabs(gq(2,1,2))**2
     &               +cdabs(gq(2,2,1))**2+cdabs(gq(2,2,2))**2)
      agsq=aveqg*fac*(cdabs(ag(1,1,1))**2+cdabs(ag(1,1,2))**2
     &               +cdabs(ag(1,2,1))**2+cdabs(ag(1,2,2))**2
     &               +cdabs(ag(2,1,1))**2+cdabs(ag(2,1,2))**2
     &               +cdabs(ag(2,2,1))**2+cdabs(ag(2,2,2))**2)

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
          if ((j > 0) .and. (k < 0)) then
            msq(j,k)=Vsq(j,k)*qasq
          elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=Vsum(k)*gasq
          elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=Vsum(k)*gqsq
          elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=Vsq(j,k)*aqsq
          elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*qgsq
         elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*agsq
          endif
      enddo
      enddo
      return
      end


