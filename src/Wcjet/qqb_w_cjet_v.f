!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_w_cjet_v(p,msq)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis
c     May, 2007.                                                  *
c***********************************************************************
c---- One-loop matrix element for W+c production, including c mass
c----averaged over initial colours and spins
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + cbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ c(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'nwz.f'
      include 'qcdcouple.f'
      include 'scheme.f'
      include 'ckm.f'
      include 'nflav.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),wprop,fac
      real(dp):: virtqg,virtgq,virtqbg,virtgqb
      real(dp):: twotDg,q(mxpart,4),dot
      complex(dp):: amp0(2,2),spp,spm,smp,smm,
     & virt_pp,virt_pm,virt_mp,virt_mm

      scheme='dred'

      msq(:,:)=0._dp

c--- calculate auxiliary momentum array - gq case
      twotDg=2._dp*dot(p,5,1)
      do k=1,4
      do j=1,5
      q(j,k)=p(j,k)
      enddo
      q(6,k)=p(5,k)-p(1,k)*mc**2/twotDg
      enddo

c---fill matrices of spinor products
      if (nwz == -1) then
         call spinoru(6,q,za,zb)
      else
         call spinoru(6,q,zb,za)
      endif

      wprop=(s(3,4)-wmass**2)**2+(wmass*wwidth)**2


      if     (nwz == -1) then
c---- basic process is g+s -> W- + c
        call tree(mc,1,2,3,4,6,amp0)
        spp=virt_pp(mc,1,2,3,4,5,q)/sqrt(wprop)
        spm=virt_pm(mc,1,2,3,4,5,q)/sqrt(wprop)
        smm=virt_mm(mc,1,2,3,4,5,q)/sqrt(wprop)
        smp=virt_mp(mc,1,2,3,4,5,q)/sqrt(wprop)

        virtgq=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &             +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

c        virtgq=virtsqwcg(2,1,3,4,5,p)/wprop   ! (with mc=mt)
      elseif (nwz == +1) then
c---- basic process is g+s~ -> W+ + c~
        call tree(mc,1,2,4,3,6,amp0)
        spp=virt_pp(mc,1,2,4,3,5,q)/sqrt(wprop)
        spm=virt_pm(mc,1,2,4,3,5,q)/sqrt(wprop)
        smm=virt_mm(mc,1,2,4,3,5,q)/sqrt(wprop)
        smp=virt_mp(mc,1,2,4,3,5,q)/sqrt(wprop)

        virtgqb=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &              +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

c        virtgqb=virtsqwcg(2,1,4,3,5,p)/wprop   ! (with mc=mt)
      else
        write(6,*) 'Problem with nwz in qqb_w_tndk.f: nwz=',nwz
        stop
      endif

c--- calculate auxiliary momentum array - qg case
      twotDg=2._dp*dot(p,5,2)
      do k=1,4
      do j=1,5
      q(j,k)=p(j,k)
      enddo
      q(6,k)=p(5,k)-p(2,k)*mc**2/twotDg
      enddo

c---fill matrices of spinor products
      if (nwz == -1) then
         call spinoru(6,q,za,zb)
      else
         call spinoru(6,q,zb,za)
      endif


      if     (nwz == -1) then
c---- basic process is s+g -> W- + c
        call tree(mc,2,1,3,4,6,amp0)
        spp=virt_pp(mc,2,1,3,4,5,q)/sqrt(wprop)
        spm=virt_pm(mc,2,1,3,4,5,q)/sqrt(wprop)
        smm=virt_mm(mc,2,1,3,4,5,q)/sqrt(wprop)
        smp=virt_mp(mc,2,1,3,4,5,q)/sqrt(wprop)

        virtqg=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &             +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

c        virtqg=virtsqwcg(1,2,3,4,5,p)/wprop   ! (with mc=mt)
      elseif (nwz == +1) then
c---- basic process is s~+g -> W+ + c~
        call tree(mc,2,1,4,3,6,amp0)
        spp=virt_pp(mc,2,1,4,3,5,q)/sqrt(wprop)
        spm=virt_pm(mc,2,1,4,3,5,q)/sqrt(wprop)
        smm=virt_mm(mc,2,1,4,3,5,q)/sqrt(wprop)
        smp=virt_mp(mc,2,1,4,3,5,q)/sqrt(wprop)

        virtqbg=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &              +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

c        virtqbg=virtsqwcg(1,2,4,3,5,p)/wprop   ! (with mc=mt)
      endif

      fac=ason2pi*gsq*gwsq**2*V*aveqg

      do j=-nflav,nflav
      do k=-nflav,nflav
      if (((j == 1) .or. (j == 3)) .and. (k == 0)) then
          msq(j,k)=fac*Vsq(j,-4)*virtqg
      elseif ((j == 0) .and. ((k == +1).or.(k == +3))) then
          msq(j,k)=fac*Vsq(-4,k)*virtgq
      elseif (((j == -1).or.(j == -3)) .and. (k == 0))then
          msq(j,k)=fac*Vsq(j,+4)*virtqbg
      elseif ((j == 0) .and. ((k == -1).or.(k == -3))) then
          msq(j,k)=fac*Vsq(+4,k)*virtgqb
      endif

      enddo
      enddo


      return
      end

