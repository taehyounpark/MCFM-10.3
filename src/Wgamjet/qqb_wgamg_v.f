!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wgamg_v(p,msq)
      implicit none
      include 'types.f'

c-----Author Keith Ellis, September 2002
c----- updated: John Campbell, August 2011 (anomalous couplings)
c----Matrix element for W gam production
c----averaged over initial colours and spins
c For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4)) + gamma(p5) + g(p6)
c For nwz=-1
c     ubar(-p1)+d(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5) + g(p6)
c---
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'zprods_decl.f'
      include 'nwz.f'
      include 'zcouple_cms.f'
      include 'scheme.f'
      include 'verbose.f'
      include 'uselfunctions.f'
      include 'blha.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: aqsq,qasq,qgsq,gqsq,agsq,gasq
      complex(dp)::qa(2,2),aq(2,2),gq(2,2),qg(2,2),ga(2,2),ag(2,2),
     & qalo(2,2),aqlo(2,2),gqlo(2,2),qglo(2,2),galo(2,2),aglo(2,2)
c     & KCD(dmax),KCC(cmax),KCB(bmax),
c     & RKCD(dmax),RKCC(cmax),RKCB(bmax),rattot,tot

      scheme='tH-V'

c flag for printing out KCheck comparison
      verbose=.false.
c flag for box functions
      uselfunctions=.false.

      call spinoru(6,p,za,zb)
c      call spinorz(6,p,za,zb)      ! For KCheck comparison

      fac=2._dp*V*gsq*abs((zesq/zxw)**2*zesq)

      if (nwz == +1) then
        call a6Wpgam_v(1,2,3,4,5,6,za,zb,qalo,qa)
        if (useblha == 0) then
        call a6Wpgam_v(2,1,3,4,5,6,za,zb,aqlo,aq)
        call a6Wpgam_v(1,6,3,4,5,2,za,zb,qglo,qg)
        call a6Wpgam_v(6,2,3,4,5,1,za,zb,galo,ga)
        call a6Wpgam_v(2,6,3,4,5,1,za,zb,gqlo,gq)
        call a6Wpgam_v(6,1,3,4,5,2,za,zb,aglo,ag)
        endif
      elseif (nwz == -1) then
        call a6Wpgam_v(2,1,4,3,5,6,zb,za,qalo,qa)
        if (useblha == 0) then
        call a6Wpgam_v(1,2,4,3,5,6,zb,za,aqlo,aq)
        call a6Wpgam_v(6,1,4,3,5,2,zb,za,qglo,qg)
        call a6Wpgam_v(2,6,4,3,5,1,zb,za,galo,ga)
        call a6Wpgam_v(6,2,4,3,5,1,zb,za,gqlo,gq)
        call a6Wpgam_v(1,6,4,3,5,2,zb,za,aglo,ag)
        endif
      endif

c--- (check of LO)
c      qasq=aveqq*fac*real(
c     & +qalo(1,1)*conjg(qalo(1,1))+qalo(1,2)*conjg(qalo(1,2))
c     & +qalo(2,1)*conjg(qalo(2,1))+qalo(2,2)*conjg(qalo(2,2)),dp)
c--- (virtual)

      qasq=aveqq*ason2pi*fac*real(
     & +qa(1,1)*conjg(qalo(1,1))+qa(1,2)*conjg(qalo(1,2))
     & +qa(2,1)*conjg(qalo(2,1))+qa(2,2)*conjg(qalo(2,2)),dp)
      aqsq=aveqq*ason2pi*fac*real(
     & +aq(1,1)*conjg(aqlo(1,1))+aq(1,2)*conjg(aqlo(1,2))
     & +aq(2,1)*conjg(aqlo(2,1))+aq(2,2)*conjg(aqlo(2,2)),dp)
      qgsq=aveqg*ason2pi*fac*real(
     & +qg(1,1)*conjg(qglo(1,1))+qg(1,2)*conjg(qglo(1,2))
     & +qg(2,1)*conjg(qglo(2,1))+qg(2,2)*conjg(qglo(2,2)),dp)
      gasq=aveqg*ason2pi*fac*real(
     & +ga(1,1)*conjg(galo(1,1))+ga(1,2)*conjg(galo(1,2))
     & +ga(2,1)*conjg(galo(2,1))+ga(2,2)*conjg(galo(2,2)),dp)
      gqsq=aveqg*ason2pi*fac*real(
     & +gq(1,1)*conjg(gqlo(1,1))+gq(1,2)*conjg(gqlo(1,2))
     & +gq(2,1)*conjg(gqlo(2,1))+gq(2,2)*conjg(gqlo(2,2)),dp)
      agsq=aveqg*ason2pi*fac*real(
     & +ag(1,1)*conjg(aglo(1,1))+ag(1,2)*conjg(aglo(1,2))
     & +ag(2,1)*conjg(aglo(2,1))+ag(2,2)*conjg(aglo(2,2)),dp)

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


