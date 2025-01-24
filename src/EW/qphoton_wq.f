!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qphoton_wq(p,msq)
      implicit none
      include 'types.f'
c----- Author John Campbell, May 2021
c----  averaged over initial colours and spins, including all crossings with a photon in initial state
c For nwz=+1
c     u(-p1)+gamma(-p2)-->W^+(n(p3)+e^+(p4)) + d(p5)
c For nwz=-1
c     ubar(-p1)+gamma(-p2)-->W^-(e^-(p3)+nbar(p4)) + dbar(p5)
c---
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ckm.f'
      include 'zprods_com.f'
      include 'nwz.f'
      include 'zcouple_cms.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),gq,qg,gqb,qbg,fac
      complex(dp):: agamtree

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(5,p,za,zb)
      fac=spinave/xn*2._dp*xn*abs((zesq/zxw)**2*zesq)

      if (nwz == -1) then
         gq=fac*(abs(agamtree(5,2,3,4,1,za,zb,-1))**2
     &          +abs(agamtree(5,2,3,4,1,za,zb,+1))**2)
         qbg=fac*(abs(agamtree(1,5,3,4,2,za,zb,-1))**2
     &           +abs(agamtree(1,5,3,4,2,za,zb,+1))**2)
         gqb=fac*(abs(agamtree(2,5,3,4,1,za,zb,-1))**2
     &           +abs(agamtree(2,5,3,4,1,za,zb,+1))**2)
         qg=fac*(abs(agamtree(5,1,3,4,2,za,zb,-1))**2
     &          +abs(agamtree(5,1,3,4,2,za,zb,+1))**2)
      elseif (nwz == +1) then
         gq=fac*(abs(agamtree(2,5,4,3,1,zb,za,-1))**2
     &          +abs(agamtree(2,5,4,3,1,zb,za,+1))**2)
         qbg=fac*(abs(agamtree(5,1,4,3,2,zb,za,-1))**2
     &           +abs(agamtree(5,1,4,3,2,zb,za,+1))**2)
         gqb=fac*(abs(agamtree(5,2,4,3,1,zb,za,-1))**2
     &           +abs(agamtree(5,2,4,3,1,zb,za,+1))**2)
         qg=fac*(abs(agamtree(1,5,4,3,2,zb,za,-1))**2
     &          +abs(agamtree(1,5,4,3,2,zb,za,+1))**2)
      endif

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize and abuse 0 to mean photon instead of gluon
      msq(j,k)=0._dp
          if     ((j == 0) .and. (k > 0)) then
            msq(j,k)=Vsum(k)*gq
          elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=Vsum(k)*gqb
          elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*qg
         elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*qbg
          endif
      enddo
      enddo

      return
      end


