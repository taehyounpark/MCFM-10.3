!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qphoton_wgamq_g(p,msq)
      implicit none
      include 'types.f'

c Gluon-photon channels only!

c-----Author Keith Ellis, May 2020
c----Matrix element for W gam production
c----averaged over initial colours and spins
c For nwz=+1
c     u(-p1)+photon(-p2)-->W^+(n(p3)+e^+(p4)) + gamma(p5) + d(p6) + g(p7)
c For nwz=-1
c     ubar(-p1)+photon(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5) + dbar(p6) + g(p7)
c---
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_com.f'
      include 'nwz.f'
      include 'zcouple_cms.f'
      include 'qcdcouple.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: ggamsq,gamgsq
      complex(dp):: gamg(2,2,2),ggam(2,2,2)

      call spinoru(7,p,za,zb)
      fac=V*gsq*abs((zesq/zxw)**2*zesq**2)

      if (nwz == -1) then
c        call a7Wpgamgamg(2,7,4,3,5,1,6,zb,za,gama)
c        call a7Wpgamgamg(7,2,4,3,5,1,6,zb,za,gamq)
c        call a7Wpgamgamg(7,1,4,3,5,2,6,zb,za,qgam)
c        call a7Wpgamgamg(1,7,4,3,5,2,6,zb,za,agam)
        call a7Wpgamgamg(6,7,4,3,5,1,2,zb,za,gamg)
        call a7Wpgamgamg(6,7,4,3,5,2,1,zb,za,ggam)
      elseif (nwz == +1) then
c        call a7Wpgamgamg(7,2,3,4,5,1,6,za,zb,gama)
c        call a7Wpgamgamg(2,7,3,4,5,1,6,za,zb,gamq)
c        call a7Wpgamgamg(1,7,3,4,5,2,6,za,zb,qgam)
c        call a7Wpgamgamg(7,1,3,4,5,2,6,za,zb,agam)
        call a7Wpgamgamg(7,6,3,4,5,1,2,za,zb,gamg) ! a(p1)+g(p2) -> W+(p34) + gam(p5) + d(p6) + u~(p7)
        call a7Wpgamgamg(7,6,3,4,5,2,1,za,zb,ggam) ! g(p1)+a(p2) -> W+(p34) + gam(p5) + d(p6) + u~(p7)
      endif

c      gamasq=spinave/xn*fac*(cdabs(gama(1,1,1))**2+cdabs(gama(1,1,2))**2
c     &                      +cdabs(gama(1,2,1))**2+cdabs(gama(1,2,2))**2
c     &                      +cdabs(gama(2,1,1))**2+cdabs(gama(2,1,2))**2
c     &                      +cdabs(gama(2,2,1))**2+cdabs(gama(2,2,2))**2)
c      gamqsq=spinave/xn*fac*(cdabs(gamq(1,1,1))**2+cdabs(gamq(1,1,2))**2
c     &                      +cdabs(gamq(1,2,1))**2+cdabs(gamq(1,2,2))**2
c     &                      +cdabs(gamq(2,1,1))**2+cdabs(gamq(2,1,2))**2
c     &                      +cdabs(gamq(2,2,1))**2+cdabs(gamq(2,2,2))**2)
c      agamsq=spinave/xn*fac*(cdabs(agam(1,1,1))**2+cdabs(agam(1,1,2))**2
c     &                      +cdabs(agam(1,2,1))**2+cdabs(agam(1,2,2))**2
c     &                      +cdabs(agam(2,1,1))**2+cdabs(agam(2,1,2))**2
c     &                      +cdabs(agam(2,2,1))**2+cdabs(agam(2,2,2))**2)
c      qgamsq=spinave/xn*fac*(cdabs(qgam(1,1,1))**2+cdabs(qgam(1,1,2))**2
c     &                      +cdabs(qgam(1,2,1))**2+cdabs(qgam(1,2,2))**2
c     &                      +cdabs(qgam(2,1,1))**2+cdabs(qgam(2,1,2))**2
c     &                      +cdabs(qgam(2,2,1))**2+cdabs(qgam(2,2,2))**2)

      ggamsq=spinave/V*fac*(cdabs(ggam(1,1,1))**2+cdabs(ggam(1,1,2))**2
     &                      +cdabs(ggam(1,2,1))**2+cdabs(ggam(1,2,2))**2
     &                      +cdabs(ggam(2,1,1))**2+cdabs(ggam(2,1,2))**2
     &                      +cdabs(ggam(2,2,1))**2+cdabs(ggam(2,2,2))**2)
      gamgsq=spinave/V*fac*(cdabs(gamg(1,1,1))**2+cdabs(gamg(1,1,2))**2
     &                      +cdabs(gamg(1,2,1))**2+cdabs(gamg(1,2,2))**2
     &                      +cdabs(gamg(2,1,1))**2+cdabs(gamg(2,1,2))**2
     &                      +cdabs(gamg(2,2,1))**2+cdabs(gamg(2,2,2))**2)

c!! Abuse of notation !!!
      msq(:,:)=0._dp
      msq(0,1) = ggamsq*2d0 ! count (u,d~) and (c,s~), for instance
      msq(1,0) = gamgsq*2d0 ! count (u,d~) and (c,s~), for instance

c      do j=-nf,nf
c      do k=-nf,nf
cc--set msq=0 to initalize
c      msq(j,k)=0._dp
c          if     ((j == 0) .and. (k < 0)) then
c!            msq(j,k)=Vsum(k)*gamasq
c          elseif ((j == 0) .and. (k > 0)) then
c!            msq(j,k)=Vsum(k)*gamqsq
c           elseif ((j > 0) .and. (k == 0)) then
c!            msq(j,k)=Vsum(j)*qgamsq
c          elseif ((j < 0) .and. (k == 0)) then
c!            msq(j,k)=Vsum(j)*agamsq
c          elseif ((j == 0) .and. (k == 0)) then
c!!!!!! Need to work out how to separate (photon, glue) from (glue, photon) .... !!!!!!
c            msq(j,k)=ggamsq*2d0 ! count (u,d~) and (c,s~), for instance
c          endif
c      enddo
c      enddo

      return
      end

      subroutine photonpdffix(fx,fxa)
      implicit none
      include 'types.f'
      real(dp):: fx(-5:5),fxa

      fx(-5:-1)=0._dp
      fx(1:5)=0._dp
      fx(1)=fxa

      return
      end

