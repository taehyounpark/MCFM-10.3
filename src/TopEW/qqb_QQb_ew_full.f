!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_QQb_ew_full(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'breit.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'first.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),ss,beta,z,
     &     qqbew(5),ggew1,ggew2,m2(-nf:nf,-nf:nf)
      integer j,k

      if(first) then
         first=.false.
         write(6,*) 'Heavy Quark mass:', mass2
      endif

c      p(1,:) = (/0._dp, 0._dp, -54.81921877743342_dp,
c     &     -54.81921877743342_dp/)
c      p(2,:) = (/0._dp, 0._dp, 5546.9288495627206_dp,
c     &     -5546.9288495627206_dp/)
c      p(3,:) = (/24.656657886284499_dp, 165.15811488069718_dp,
c     &     -5267.1711511618723_dp, 5272.6466199927481_dp/)
c      p(4,:) = (/-24.656657886284499_dp, -165.15811488069718_dp,
c     &     -224.93847962341448_dp, 329.10144834740561_dp/)

      call dotem(4,p,s)

      ss = s(1,2)
      beta = sqrt(1._dp-4._dp*mt**2/ss)
      z = (s(2,3)-s(1,3))/ss/beta

c      print*, beta, ss

      call qqbQQb_ew_oneloop(qqbew,ss,beta,z)
      call ggQQb_ew_oneloop(ggew1,ss,beta,z)
      call ggQQb_ew_oneloop(ggew2,ss,beta,-z)
      call qqb_QQb(p,m2)

      msq = 0._dp

      do j=-nf,nf
         k=-j
         if((j == 0) .and. (k == 0)) then
            msq(j,k) = ggew1 + ggew2
         else if((j > 0) .and. (k < 0)) then
            msq(j,k) = qqbew(j)
         else if((j < 0 ) .and. (k > 0)) then
            msq(j,k) = qqbew(k)
         end if
      end do

      msq = msq*m2!/gsq**2*16._dp*pi**2*0.01_dp
c      msq(0,0) = m2(0,0)/gsq**2*16._dp*pi**2*0.01_dp

c      print*, 'phase space point: '
c      print*,
c      print*, 'p(1,:) ', p(1,:)
c      print*, 'p(2,:) ', p(2,:)
c      print*, 'p(3,:) ', p(3,:)
c      print*, 'p(4,:) ', p(4,:)
c      print*, 'beta: ', beta
c      print*, 's: ', ss
c      print*, 'z: ', z
c      print*, "dd~ :", msq(1,-1)
c      print*, "dd~(born) :", m2(1,-1)/gsq**2*16._dp*pi**2*0.01_dp
c      print*, "uu~ :", msq(2,-2)
c      print*, "uu~(born) :", m2(2,-2)/gsq**2*16._dp*pi**2*0.01_dp
c      print*, "gg(born) :", m2(0,0)
c      stop

      end subroutine qqb_QQb_ew_full
