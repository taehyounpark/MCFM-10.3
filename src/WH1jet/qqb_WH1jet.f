!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_WH1jet(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c---for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+b(p6)
c---for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> b(p5)+b(p6)
c   for the moment --- radiation only from initial line

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'hbbparams.f'
      include 'hdecaymode.f'
      include 'blha.f'
      integer:: j,k
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: radi,radi_ww,hdecay
      real(dp):: qqbWHg,qbqWHg,qgWHq,gqWHq,gqbWHqb,qbgWHqb

      msq(:,:)=zip

      if (hdecaymode == 'wpwm') then
        call dotem(9,p,s)
        call hwwdecay(p,5,6,7,8,hdecay)
        qqbWHg=aveqq*radi_ww(1,2,9,5,6,7,8,3,4)*hdecay
        qbqWHg=aveqq*radi_ww(2,1,9,5,6,7,8,3,4)*hdecay
        qgWHq=-radi_ww(1,9,2,5,6,7,8,3,4)*aveqg*hdecay
        gqWHq=-radi_ww(2,9,1,5,6,7,8,3,4)*aveqg*hdecay
        gqbWHqb=-radi_ww(9,2,1,5,6,7,8,3,4)*aveqg*hdecay
        qbgWHqb=-radi_ww(9,1,2,5,6,7,8,3,4)*aveqg*hdecay
      else
        call dotem(7,p,s)
c selection of other Higgs decay modes in function radi
        if ((useblha == 0).or.((blhafl(1) /= 0).and.(blhafl(2) /= 0))) then
        qqbWHg=aveqq*radi(1,2,7,5,6,3,4)
        endif
        if (useblha == 0) then
        qbqWHg=aveqq*radi(2,1,7,5,6,3,4)
        endif
        if ((useblha == 0).or.((blhafl(1) /= 0).and.(blhafl(2) == 0))) then
        qgWHq=-radi(1,7,2,5,6,3,4)*aveqg
        endif
        if (useblha == 0) then
        gqWHq=-radi(2,7,1,5,6,3,4)*aveqg
        endif
        if ((useblha == 0).or.((blhafl(1) == 0).and.(blhafl(2) /= 0))) then
        gqbWHqb=-radi(7,2,1,5,6,3,4)*aveqg
        endif
        if (useblha == 0) then
        qbgWHqb=-radi(7,1,2,5,6,3,4)*aveqg
        endif
      endif
      if (p(7,4)<0.) then
         qqbWHg=-qqbWHg
      endif

c      write(6,*) 'qqbWHg',qqbWHg
c      write(6,*) 'qbqWHg',qbqWHg
c      write(6,*) 'qbgWHqb',qbgWHqb
c      write(6,*) 'gqbWHqb',gqbWHqb
c      write(6,*) 'qgWHq',qgWHq
c      write(6,*) 'gqWHq',gqWHq

      do j=-nf,nf
      do k=-nf,nf

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*qqbWHg
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*qbqWHg
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWHq
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWHqb
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWHq
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWHqb
      endif

      enddo
      enddo

c--- adjust for fixed H->bb BR if necessary
      if ((hdecaymode == 'bqba') .and. (FixBrHbb)) then
        msq(:,:)=msq(:,:)*GamHbb/GamHbb0
      endif

      return
      end


