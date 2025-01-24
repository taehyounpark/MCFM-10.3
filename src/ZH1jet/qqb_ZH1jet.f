!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_ZH1jet(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z +g(p7)
c                           |    |
c                           |    --> fermion(p3)+antifermion(p4)
c                           |
c                           ---> b(p5)+b(p6)
c   for the moment --- radiation only from initial line
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'hbbparams.f'
      include 'hdecaymode.f'
      include 'zcouple_cms.f'
      include 'blha.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: radiLL,radiLL_ww
      real(dp):: qqbZHgL,qbqZHgL,qgZHqL,gqZHqL,gqbZHqbL,qbgZHqbL
      real(dp):: qqbZHgR,qbqZHgR,qgZHqR,gqZHqR,gqbZHqbR,qbgZHqbR


      msq(:,:)=zip

      if (hdecaymode=='wpwm') then
         call dotem(9,p,s)
         qqbZHgL=aveqq*radiLL_ww(1,2,9,5,6,7,8,3,4)
         qqbZHgR=aveqq*radiLL_ww(1,2,9,5,6,7,8,4,3)
         qbqZHgL=aveqq*radiLL_ww(2,1,9,5,6,7,8,3,4)
         qbqZHgR=aveqq*radiLL_ww(2,1,9,5,6,7,8,4,3)

         qgZHqL=-radiLL_ww(1,9,2,5,6,7,8,3,4)*aveqg
         qgZHqR=-radiLL_ww(1,9,2,5,6,7,8,4,3)*aveqg
         gqZHqL=-radiLL_ww(2,9,1,5,6,7,8,3,4)*aveqg
         gqZHqR=-radiLL_ww(2,9,1,5,6,7,8,4,3)*aveqg

         gqbZHqbL=-radiLL_ww(9,2,1,5,6,7,8,3,4)*aveqg
         gqbZHqbR=-radiLL_ww(9,2,1,5,6,7,8,4,3)*aveqg

         qbgZHqbL=-radiLL_ww(9,1,2,5,6,7,8,3,4)*aveqg
         qbgZHqbR=-radiLL_ww(9,1,2,5,6,7,8,4,3)*aveqg
      else
         call dotem(7,p,s)

         if ((useblha == 0).or.((blhafl(1) /= 0).and.(blhafl(2) /= 0))) then
         qqbZHgL=aveqq*radiLL(1,2,7,5,6,3,4)
         qqbZHgR=aveqq*radiLL(1,2,7,5,6,4,3)
         endif
         if (useblha == 0) then
         qbqZHgL=aveqq*radiLL(2,1,7,5,6,3,4)
         qbqZHgR=aveqq*radiLL(2,1,7,5,6,4,3)
         endif

         if ((useblha == 0).or.((blhafl(1) /= 0).and.(blhafl(2) == 0))) then
         qgZHqL=-radiLL(1,7,2,5,6,3,4)*aveqg
         qgZHqR=-radiLL(1,7,2,5,6,4,3)*aveqg
         endif
         if (useblha == 0) then
         gqZHqL=-radiLL(2,7,1,5,6,3,4)*aveqg
         gqZHqR=-radiLL(2,7,1,5,6,4,3)*aveqg
         endif

         if ((useblha == 0).or.((blhafl(1) == 0).and.(blhafl(2) /= 0))) then
         gqbZHqbL=-radiLL(7,2,1,5,6,3,4)*aveqg
         gqbZHqbR=-radiLL(7,2,1,5,6,4,3)*aveqg
         endif

         if (useblha == 0) then
         qbgZHqbL=-radiLL(7,1,2,5,6,3,4)*aveqg
         qbgZHqbR=-radiLL(7,1,2,5,6,4,3)*aveqg
         endif
      endif
      if (p(7,4)<0.) then
         qqbZHgL=-qqbZHgL
         qqbZHgR=-qqbZHgR
      endif

      do j=-nf,nf
      do k=-nf,nf

      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) cycle
      if ((j == 0) .and. (k == 0)) cycle

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=(abs(zL(j)*zl1)**2+abs(zR(j)*zr1)**2)*qqbZHgL
     &            +(abs(zR(j)*zl1)**2+abs(zL(j)*zr1)**2)*qqbZHgR
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=(abs(zL(k)*zl1)**2+abs(zR(k)*zr1)**2)*qbqZHgL
     &            +(abs(zR(k)*zl1)**2+abs(zL(k)*zr1)**2)*qbqZHgR
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=(abs(zL(j)*zl1)**2+abs(zR(j)*zr1)**2)*qgZHqL
     &            +(abs(zR(j)*zl1)**2+abs(zL(j)*zr1)**2)*qgZHqR
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=(abs(zL(-j)*zl1)**2+abs(zR(-j)*zr1)**2)*qbgZHqbL
     &            +(abs(zR(-j)*zl1)**2+abs(zL(-j)*zr1)**2)*qbgZHqbR
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=(abs(zL(k)*zl1)**2+abs(zR(k)*zr1)**2)*gqZHqL
     &            +(abs(zR(k)*zl1)**2+abs(zL(k)*zr1)**2)*gqZHqR
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=(abs(zL(-k)*zl1)**2+abs(zR(-k)*zr1)**2)*gqbZHqbL
     &            +(abs(zR(-k)*zl1)**2+abs(zL(-k)*zr1)**2)*gqbZHqbR
      endif

      enddo
      enddo
c---  adjust for fixed H->bb BR if necessary
      if ((hdecaymode == 'bqba') .and. (FixBrHbb)) then
         msq(:,:)=msq(:,:)*GamHbb/GamHbb0
      endif


      return
      end


