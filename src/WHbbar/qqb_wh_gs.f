!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wh_gs(p,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c---for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+bb(p6)
c---for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> b(p5)+bb(p6)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer:: j,k,nd

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & sub17_2(4),sub27_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      external qqb_wh,donothing_gvec

      ndmax=2

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     & qqb_wh,donothing_gvec)
      call dips(2,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     & qqb_wh,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo

      if  ((j == 0) .and. (k == 0)) then
         goto 20
      elseif  ((j > 0) .and. (k < 0)
     &     .or.(j < 0) .and. (k > 0)) then
         msq(1,j,k)=2._dp*cf*sub17_2(qq)*msq17_2(j,k)
         msq(2,j,k)=2._dp*cf*sub27_1(qq)*msq27_1(j,k)

      elseif ((j  /=  0) .and. (k == 0)) then
         msq(1,j,k)=0._dp
         msq(2,j,k)=2._dp*tr*sub27_1(qg)*(msq27_1(j,+1)+msq27_1(j,+2)
     &   +msq27_1(j,+3)+msq27_1(j,+4)+msq27_1(j,+5)
     &                                 +msq27_1(j,-1)+msq27_1(j,-2)
     &   +msq27_1(j,-3)+msq27_1(j,-4)+msq27_1(j,-5))

      elseif ((j == 0) .and. (k  /=  0)) then
         msq(1,j,k)=2._dp*tr*sub17_2(qg)*(msq17_2(+1,k)+msq17_2(+2,k)
     &   +msq17_2(+3,k)+msq17_2(+4,k)+msq17_2(+5,k)
     &                                 +msq17_2(-1,k)+msq17_2(-2,k)
     &   +msq17_2(-3,k)+msq17_2(-4,k)+msq17_2(-5,k))
      msq(2,j,k)=0._dp

      endif
 20   continue
      enddo
      enddo

      return
      end


